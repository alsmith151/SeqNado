"""
User input module for SeqNado configuration using the new Pydantic models.
"""

import json
import os
import pathlib
import sys
from datetime import datetime
from typing import Dict, List, Optional, Any
from pathlib import Path
from loguru import logger
from pydantic import ValidationError
import jinja2

from seqnado import (
    Assay,
    SpikeInMethod,
    PileupMethod,
    PeakCallingMethod,
    QuantificationMethod,
    SNPCallingMethod,
    MethylationMethod
)
from seqnado.config import (
    SeqnadoConfig,
    GenomeConfig,
    ProjectConfig,
    BigwigConfig,
    PlottingConfig,
    PeakCallingConfig,
    SpikeInConfig,
    UCSCHubConfig,
    RNAQuantificationConfig,
    SNPCallingConfig,
    MCCConfig,
    MethylationConfig,
    MLDatasetConfig,
    UserFriendlyError,
    ATACAssayConfig,
    ChIPAssayConfig,
    CATAssayConfig,
    RNAAssayConfig,
    SNPAssayConfig,
    MCCAssayConfig,
    MethylationAssayConfig,
    CRISPRAssayConfig,
    BowtieIndex,
    STARIndex,
    AssaySpecificConfig,
)


def get_user_input(
    prompt: str,
    default: Optional[str] = None,
    is_boolean: bool = False,
    choices: Optional[List[str]] = None,
    is_path: bool = False,
    required: bool = True,
) -> str:
    """
    Prompt the user for input with validation for choices, boolean values, or path existence.
    Re-prompts until valid input is provided.
    """
    while True:
        # Construct the prompt suffix based on choices or default
        if choices:
            prompt_suffix = f"({'/'.join(choices)})"
        elif default is not None:
            prompt_suffix = f"(default: {default})"
        else:
            prompt_suffix = ""

        user_input = input(f"{prompt} {prompt_suffix}: ").strip()

        # Handle empty input and apply default if available
        if not user_input:
            if default is not None:
                user_input = default
            elif not required:
                return None
            else:
                print("Input cannot be empty. Please try again.")
                continue

        # Validate boolean input
        if is_boolean:
            if user_input.lower() in {"yes", "y", "true", "1"}:
                return True
            elif user_input.lower() in {"no", "n", "false", "0"}:
                return False
            else:
                print(
                    "Invalid boolean value. Please enter yes/no, y/n, true/false, or 1/0."
                )
                continue

        # Validate against allowed choices
        if choices and user_input not in choices:
            print(f"Invalid choice. Please choose from: {', '.join(choices)}")
            continue

        # Validate path existence if required
        if is_path and user_input and not os.path.exists(user_input):
            print(f"The path '{user_input}' does not exist. Please try again.")
            continue

        return user_input


def load_genome_configs(assay: Assay) -> Dict[str, BowtieIndex | STARIndex]:
    """Load genome configurations from the config file."""
    config_path = (
        pathlib.Path(os.getenv("SEQNADO_CONFIG", pathlib.Path.home()))
        / ".config/seqnado/genome_config.json"
    )
    if not config_path.exists():
        logger.error("Genome config not found. Run 'seqnado-init' first.")
        sys.exit(1)

    with open(config_path) as f:
        genome_data = json.load(f)

    # Convert to GenomeConfig objects
    genome_configs = {}
    for name, config_dict in genome_data.items():
        genome_index = BowtieIndex if assay != Assay.RNA else STARIndex
        index_key = "bt2_index" if assay != Assay.RNA else "star_index"
        if index_key not in config_dict:
            print(
                f"\n❌ Configuration error for genome '{name}': Missing '{index_key}' key.\n"
            )
            sys.exit(1)

        try:
            genome_configs[name] = GenomeConfig(
                name=name,
                index=genome_index(prefix=config_dict.get(index_key)),
                fasta=pathlib.Path(config_dict["fasta"]),
                chromosome_sizes=pathlib.Path(config_dict["chromosome_sizes"])
                if config_dict.get("chromosome_sizes")
                else None,
                gtf=pathlib.Path(config_dict["gtf"])
                if config_dict.get("gtf")
                else None,
                genes=pathlib.Path(config_dict["genes"])
                if config_dict.get("genes")
                else None,
                blacklist=pathlib.Path(config_dict["blacklist"])
                if config_dict.get("blacklist")
                else None,
                organism=config_dict.get("organism"),
                version=config_dict.get("version"),
            )
        except ValidationError as e:
            for error in e.errors():
                if error.get("type") == "value_error" and hasattr(
                    error.get("ctx", {}), "get"
                ):
                    # Extract the custom error message if it exists
                    error_msg = str(error.get("ctx", {}).get("error", ""))
                    if "❌" in error_msg:  # Our custom formatted error
                        print(f"\n{error_msg}\n")
                        sys.exit(1)

            # Fallback for other validation errors
            print(f"\n❌ Configuration error for genome '{name}':")
            print(f"   {str(e)}\n")
            sys.exit(1)
        except UserFriendlyError as e:
            # Display our clean error message without traceback
            print(f"\n{str(e)}\n")
            sys.exit(1)

    return genome_configs


def get_project_config() -> ProjectConfig:
    """Get project configuration from user input."""
    username = os.getenv("USER", "unknown_user")
    today = datetime.strftime(datetime.today(), "%Y-%m-%d")

    project_name = get_user_input(
        "Project name?", default=f"{username}_project"
    ).replace(" ", "_")

    project_dir = get_user_input(
        "Project directory?", default=os.getcwd(), is_path=False
    )

    return ProjectConfig(
        name=project_name, date=today, directory=pathlib.Path(project_dir)
    )


def select_genome_config(genome_configs: Dict[str, GenomeConfig]) -> GenomeConfig:
    """Allow user to select a genome configuration."""
    available_genomes = list(genome_configs.keys())

    genome_name = get_user_input(
        f"Genome? (Available: {', '.join(available_genomes)})",
        default=available_genomes[0] if available_genomes else "",
    )

    while genome_name not in genome_configs:
        print(
            f"Genome '{genome_name}' is not configured. Please choose from: {', '.join(available_genomes)}"
        )
        genome_name = get_user_input(
            "Genome?", default=available_genomes[0] if available_genomes else ""
        )

    return genome_configs[genome_name]


def get_bigwig_config() -> Optional[BigwigConfig]:
    """Get bigwig configuration if user wants to create bigwigs."""
    make_bigwigs = get_user_input("Make Bigwigs?", default="no", is_boolean=True)

    if not make_bigwigs:
        return None

    pileup_method = get_user_input(
        "Bigwig method:", choices=[m.value for m in PileupMethod], default="deeptools"
    )

    binsize = get_user_input("Binsize for bigwigs:", default="10", required=False)

    return BigwigConfig(pileup_method=[PileupMethod(pileup_method)], binsize=binsize)


def get_plotting_config() -> Optional[PlottingConfig]:
    """Get plotting configuration if user wants plotting."""
    perform_plotting = get_user_input(
        "Perform plotting?", default="no", is_boolean=True
    )

    if not perform_plotting:
        return None

    coordinates = get_user_input(
        "Path to bed file with coordinates for plotting:", required=False, is_path=True
    )

    genes = get_user_input("Path to bed file with genes:", required=False, is_path=True)

    return PlottingConfig(coordinates=coordinates, genes=genes)


def get_peak_calling_config(assay: Assay) -> Optional[PeakCallingConfig]:
    """Get peak calling configuration for assays that support it."""
    if assay not in [Assay.CHIP, Assay.ATAC, Assay.CAT]:
        return None

    call_peaks = get_user_input("Call peaks?", default="yes", is_boolean=True)

    if not call_peaks:
        return None

    # Set default based on assay
    default_method = {
        Assay.CHIP: "lanceotron",
        Assay.ATAC: "lanceotron",
        Assay.CAT: "seacr",
    }.get(assay, "lanceotron")

    peak_calling_method = get_user_input(
        "Peak calling method:",
        choices=[m.value for m in PeakCallingMethod],
        default=default_method,
    )

    consensus_counts = get_user_input(
        "Generate consensus counts from Design merge column?",
        default="no",
        is_boolean=True,
    )

    return PeakCallingConfig(
        method=[PeakCallingMethod(peak_calling_method)],
        consensus_counts=consensus_counts,
    )


def get_spikein_config(assay: Assay) -> Optional[SpikeInConfig]:
    """Get spike-in configuration for assays that support it."""
    if assay not in [Assay.CHIP, Assay.CAT]:
        return None

    spikein = get_user_input("Do you have spikein?", default="no", is_boolean=True)

    if not spikein:
        return None

    normalisation_method = get_user_input(
        "Normalisation method:",
        choices=[m.value for m in SpikeInMethod],
        default="orlando",
    )

    reference_genome = get_user_input("Reference genome:", default="hg38")
    spikein_genome = get_user_input("Spikein genome:", default="dm6")

    return SpikeInConfig(
        method=SpikeInMethod(normalisation_method),
        endogenous_genome=reference_genome,
        exogenous_genome=spikein_genome,
    )


def get_ucsc_hub_config() -> Optional[UCSCHubConfig]:
    """Get UCSC hub configuration if user wants it."""
    make_ucsc_hub = get_user_input("Make UCSC hub?", default="no", is_boolean=True)

    if not make_ucsc_hub:
        return None

    username = os.getenv("USER", "unknown_user")

    directory = get_user_input("UCSC hub directory:", default="seqnado_output/hub/")

    email = get_user_input(
        "What is your email address?", default=f"{username}@example.com"
    )

    color_by = get_user_input("Color by (for UCSC hub):", default="samplename")

    return UCSCHubConfig(directory=directory, email=email, color_by=color_by)


def get_ml_dataset_config(assay: Assay) -> Optional[MLDatasetConfig]:
    """Get ML dataset configuration for supported assays."""
    if assay not in [Assay.ATAC, Assay.CHIP, Assay.CAT]:
        return None

    make_dataset = get_user_input("Make dataset for ML?", default="no", is_boolean=True)

    if not make_dataset:
        return None

    use_regions = get_user_input(
        "Use regions BED file?", default="yes", is_boolean=True
    )

    if use_regions:
        regions_bed = get_user_input(
            "Path to regions BED file:", default="path/to/regions.bed", is_path=True
        )
        return MLDatasetConfig(regions_bed=pathlib.Path(regions_bed))
    else:
        binsize = int(get_user_input("Binsize for dataset:", default="1000"))
        return MLDatasetConfig(binsize=binsize)


def get_rna_quantification_config() -> Optional[RNAQuantificationConfig]:
    """Get RNA quantification configuration."""
    method = get_user_input(
        "Quantification method:",
        choices=[m.value for m in QuantificationMethod],
        default="feature_counts",
    )

    salmon_index = None
    if method == "salmon":
        salmon_index = get_user_input(
            "Salmon index path:", default="path/to/salmon_index"
        )

    run_deseq2 = get_user_input("Run DESeq2?", default="no", is_boolean=True)

    return RNAQuantificationConfig(
        method=QuantificationMethod(method),
        salmon_index=salmon_index,
        run_deseq2=run_deseq2,
    )


def get_snp_calling_config() -> Optional[SNPCallingConfig]:
    """Get SNP calling configuration."""
    call_snps = get_user_input("Call SNPs?", default="no", is_boolean=True)

    if not call_snps:
        return None

    method = get_user_input(
        "SNP calling method:",
        choices=[m.value for m in SNPCallingMethod],
        default="bcftools",
    )

    annotate_snps = get_user_input("Annotate SNPs?", default="no", is_boolean=True)

    snp_database = None
    if annotate_snps:
        snp_database = get_user_input(
            "Path to SNP database:", default="path/to/snp_database", is_path=True
        )

    return SNPCallingConfig(
        method=SNPCallingMethod(method),
        annotate_snps=annotate_snps,
        snp_database=snp_database,
    )


def get_mcc_config() -> Optional[MCCConfig]:
    """Get MCC (Capture-C) configuration."""
    viewpoints = get_user_input(
        "Path to viewpoints file:", default="path/to/viewpoints.bed", is_path=True
    )

    resolutions_str = get_user_input(
        "Resolutions for MCC cooler files (comma-separated):", default="100,1000"
    )

    resolutions = [int(r.strip()) for r in resolutions_str.split(",")]

    return MCCConfig(viewpoints=pathlib.Path(viewpoints), resolutions=resolutions)


def get_methylation_config() -> Optional[MethylationConfig]:
    """Get methylation calling configuration."""
    call_methylation = get_user_input(
        "Call methylation?", default="no", is_boolean=True
    )

    if not call_methylation:
        return None

    methylation_assay = get_user_input(
        "Methylation assay:",
        choices=[m.value for m in MethylationMethod],
        default="taps",
    )

    return MethylationConfig(assay=MethylationMethod(methylation_assay))


def build_assay_config(
    assay: Assay, genome_config: GenomeConfig
) -> Optional[AssaySpecificConfig]:
    """Build assay-specific configuration based on the assay type."""

    # Get common configurations
    bigwigs = get_bigwig_config()
    plotting = get_plotting_config()
    ucsc_hub = get_ucsc_hub_config()
    create_heatmaps = get_user_input("Make heatmaps?", default="no", is_boolean=True)
    geo_files = get_user_input(
        "Generate GEO submission files?", default="no", is_boolean=True
    )

    base_config = {
        "genome": genome_config,
        "bigwigs": bigwigs,
        "plotting": plotting,
        "ucsc_hub": ucsc_hub,
        "create_heatmaps": create_heatmaps,
        "create_geo_submission_files": geo_files,
    }

    
    match assay:
        case Assay.ATAC:
            tn5_shift = get_user_input("Shift ATAC reads?", default="yes", is_boolean=True)
            peak_calling = get_peak_calling_config(assay)
            dataset_for_ml = get_ml_dataset_config(assay)

            return ATACAssayConfig(
                **base_config,
                tn5_shift=tn5_shift,
                peak_calling=peak_calling,
                dataset_for_ml=dataset_for_ml,
            )

        case Assay.CHIP:
            spikein = get_spikein_config(assay)
            peak_calling = get_peak_calling_config(assay)
            dataset_for_ml = get_ml_dataset_config(assay)

            return ChIPAssayConfig(
                **base_config,
                spikein=spikein,
                peak_calling=peak_calling,
                dataset_for_ml=dataset_for_ml,
            )

        case Assay.CAT:
            tn5_shift = get_user_input("Shift CAT reads?", default="no", is_boolean=True)
            spikein = get_spikein_config(assay)
            peak_calling = get_peak_calling_config(assay)
            dataset_for_ml = get_ml_dataset_config(assay)

            return CATAssayConfig(
                **base_config,
                tn5_shift=tn5_shift,
                spikein=spikein,
                peak_calling=peak_calling,
                dataset_for_ml=dataset_for_ml,
            )

        case Assay.RNA:
            rna_quantification = get_rna_quantification_config()

            return RNAAssayConfig(**base_config, rna_quantification=rna_quantification)

        case Assay.SNP:
            snp_calling = get_snp_calling_config()

            return SNPAssayConfig(**base_config, snp_calling=snp_calling)

        case Assay.MCC:
            mcc = get_mcc_config()
            return MCCAssayConfig(**base_config, mcc=mcc)

        case Assay.METH:
            methylation = get_methylation_config()

            return MethylationAssayConfig(**base_config, methylation=methylation)

        case Assay.CRISPR:
            return CRISPRAssayConfig(**base_config)

        case _:
            raise ValueError(f"Unsupported assay type: {assay}")


def build_workflow_config(assay: Assay, seqnado_version: str) -> SeqnadoConfig:
    """Build complete workflow configuration from user input."""

    logger.info(
        f"Building workflow configuration for {assay.value} assay with SeqNado version {seqnado_version}"
    )
    logger.info(f"Current Working Directory: {os.getcwd()}")
    logger.debug(f"Assay: {assay}, SeqNado version: {seqnado_version}")

    # Load available genome configurations
    genome_configs = load_genome_configs(assay)

    # Get project configuration
    project = get_project_config()

    # Select genome configuration
    genome = select_genome_config(genome_configs)

    # Get metadata path
    metadata_path = get_user_input(
        "Path to metadata file:", default="metadata.csv", is_path=True
    )

    # Build assay-specific configuration
    assay_config = build_assay_config(assay, genome)

    try:
        workflow_config = SeqnadoConfig(
            assay=assay,
            project=project,
            genome=genome,
            metadata=pathlib.Path(metadata_path),
            assay_config=assay_config,
        )
        return workflow_config

    except ValidationError as e:
        logger.error(f"Configuration validation error: {e}")
        sys.exit(1)


def render_config(
    template: Path,
    workflow_config: SeqnadoConfig,
    outfile: Path,
    all_options: bool = False,
) -> None:
    """Render the workflow configuration to a file."""

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(template.parent))
    template = env.get_template(template.name)

    # Convert the Pydantic model to a dictionary for rendering
    config_dict = workflow_config.model_dump(mode="json", exclude_none=not all_options)

    try:
        rendered_content = template.render(**config_dict)
    except jinja2.TemplateError as e:
        logger.error(f"Template rendering error: {e}")
        sys.exit(1)

    with open(outfile, "w") as f:
        f.write(rendered_content)
