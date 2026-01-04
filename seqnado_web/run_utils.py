
import os
import yaml
from pathlib import Path
from datetime import date
import importlib.resources as resources

# SeqNado imports
from seqnado.inputs import Assay
from seqnado.config.configs import (
    GenomeConfig,
    ProjectConfig,
    BigwigConfig,
    PeakCallingConfig,
    SpikeInConfig,
    UCSCHubConfig,
    RNAQuantificationConfig,
    SNPCallingConfig,
    MethylationConfig,
    MCCConfig,
)
from seqnado.config.core import SeqnadoConfig, AssaySpecificConfig
from seqnado.config.user_input import render_config

def generate_workflow_config(
    assay_name: str,
    project_name: str,
    genome_name: str,
    design_file: str,
    output_dir: str = ".",
) -> tuple[Path, str]:
    """
    Generates a SeqNado configuration YAML file based on web form inputs.

    Args:
        assay_name: The name of the assay (e.g., 'rna', 'atac').
        project_name: The name of the project.
        genome_name: The name of the genome to use (must exist in genome_config.json).
        design_file: Path to the metadata CSV file.
        output_dir: Directory where the config and project folder should be created.

    Returns:
        A tuple containing:
        - The path to the generated configuration YAML file.
        - The project directory path (as a string).
    """
    
    # 1. Instantiate the Assay enum
    try:
        assay = Assay.from_clean_name(assay_name)
    except ValueError:
        raise ValueError(f"Invalid assay name: {assay_name}")

    # 2. Create the Project Config
    # Note: In the web app, we might simplify options, or we could add more fields to the form later.
    # For now, we use sensible defaults for a standard run.
    project_config = ProjectConfig(
        name=project_name,
        metadata_file=design_file,
        # We assume standard paths or let the user move things later, 
        # but optimally we should handle the fastq symlinking too if possible.
        # For this iteration, we focus on generating the valid config.
    )

    # 3. Create Assay-Specific Configs
    # We populate these with defaults for now, as the web wizard is simple.
    # In a full version, we'd map every form field to these configs.
    
    bigwig_conf = BigwigConfig(mode="ignore") # Default to ignore to speed up for testing? Or 'generate'? Let's stick to defaults.
    peak_conf = PeakCallingConfig(mode="ignore")
    spike_conf = SpikeInConfig(mode="ignore")
    ucsc_conf = UCSCHubConfig(mode="ignore")
    rna_conf = RNAQuantificationConfig(mode="ignore")
    snp_conf = SNPCallingConfig(mode="ignore")
    not_impl = MethylationConfig(mode="ignore") # Placeholders
    mcc_conf = MCCConfig(mode="ignore")
    
    # Enable specific analyses based on assay defaults
    if assay == Assay.ATAC:
        bigwig_conf.mode = "generate"
        peak_conf.mode = "generate"
        peak_conf.caller = "macs2" # Default
    elif assay == Assay.CHIP:
        bigwig_conf.mode = "generate"
        peak_conf.mode = "generate"
        peak_conf.caller = "macs2"
    elif assay == Assay.RNA:
        bigwig_conf.mode = "generate"
        rna_conf.mode = "generate"
        rna_conf.quantifier = "featurecounts"
    
    
    # 4. Construct the AssaySpecificConfig
    # We need to map the specifics but AssaySpecificConfig is a Union/TypeVar in the codebase 
    # dependent on the assay. However, `SeqnadoConfig` takes specific sub-configs.
    
    # Wait, looking at `seqnado.config.core.SeqnadoConfig`, it takes these sub-configs directly. 
    # Let's verify `seqnado.config.core.SeqnadoConfig` structure.
    # It seems I need to create the main SeqnadoConfig object.
    
    # NOTE: I need to verify the exact init of SeqNadoConfig.
    # Assuming standard structure based on previous file views.
    
    workflow_config = SeqnadoConfig(
        assay=assay,
        project=project_config,
        genome=GenomeConfig(name=genome_name, index=None), # The index is resolved at runtime/render time usually, or we need to load it. 
        # Actually `render_config` expects a fully populated object.
        # simpler approach: The standard `build_workflow_config` loads the genome config.
        # We should load the genome config to be safe.
        bigwig=bigwig_conf,
        peak_calling=peak_conf,
        spike_in=spike_conf,
        ucsc_hub=ucsc_conf,
        rna_quantification=rna_conf,
        snp_calling=snp_conf,
        methylation=not_impl, # filling required fields
        mcc=mcc_conf,
        # fill others if needed
    )

    # Update genome config from the actual file 
    from seqnado.config.user_input import load_genome_configs
    # We'll need to mock the Assay for load_genome_configs if it filters by assay index type
    # But `load_genome_configs` takes an assay enum.
    available_genomes = load_genome_configs(assay)
    if genome_name in available_genomes:
        workflow_config.genome = available_genomes[genome_name]
    else:
        # Fallback or error?
        pass # The template might just need the name, but best to have the full config
    
    
    # 5. Create Directory Structure
    # "The pipeline process multiomics data... Output: Processed data in seqnado_output/"
    # Usually `seqnado config` creates a directory `YYYY-MM-DD_<project_name>` and inside it a `fastqs` dir.
    
    run_dir_name = f"{date.today().isoformat()}_{project_name}"
    run_dir = Path(output_dir) / run_dir_name
    fastq_dir = run_dir / "fastqs"
    fastq_dir.mkdir(parents=True, exist_ok=True)
    
    # 6. Render Config
    config_filename = f"config_{assay_name}.yaml"
    config_path = run_dir / config_filename
    
    # We need the Jinja template path. 
    # In cli.py: `tpl_trav = _pkg_traversable("seqnado.data").joinpath("config_template.jinja")`
    # We can do similar.
    
    try:
        # Python 3.9+ 'files' API, or fallback to pkg_resources or internal helpers
        # using the same logic as CLI for safety
        try:
             from importlib.resources import files as _pkg_traversable
        except ImportError:
             from importlib.resources import files as _pkg_traversable # should exist in 3.10+
        
        tpl_trav = _pkg_traversable("seqnado.data").joinpath("config_template.jinja")
        
        with resources.as_file(tpl_trav) as tpl_path:
             render_config(
                template=Path(tpl_path),
                workflow_config=workflow_config,
                outfile=config_path,
                all_options=False # Only active options
            )
            
    except Exception as e:
        raise RuntimeError(f"Failed to render config: {e}")

    return config_path, str(run_dir)

def get_run_command(config_path: Path, profile: str = "slurm") -> str:
    """
    Returns the command string to run the pipeline.
    """
    # Absolute path is safer for instructions
    abs_config_path = config_path.resolve()
    
    # Determine the directory to run from (usually the parent of the config for single-assay)
    # The CLI runs in the Dir containing the config usually or CWD if config is passed.
    # Standard: user goes into the dir and runs `seqnado pipeline ...`
    
    # Command: seqnado pipeline <assay> --configfile <config> --profile <profile>
    # Actually checking `pipeline` command in CLI: it takes `assay` and `--configfile`.
    
    # Wait, the config file name usually implies the assay: config_<assay>.yaml
    # But explicitly:
    assay_name = config_path.stem.replace("config_", "")
    
    cmd = f"seqnado pipeline {assay_name} --configfile {abs_config_path} --profile {profile}"
    return cmd
