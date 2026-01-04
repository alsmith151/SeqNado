import os
import shutil
from pathlib import Path
from datetime import date
import importlib.resources as resources
from typing import List, Dict

from seqnado.inputs import Assay
from seqnado.config.user_input import render_config
from seqnado.config.core import SeqnadoConfig, AssaySpecificConfig
from seqnado.config.configs import (
    GenomeConfig, ProjectConfig, BigwigConfig, PeakCallingConfig, 
    SpikeInConfig, UCSCHubConfig, RNAQuantificationConfig, 
    SNPCallingConfig, MethylationConfig, MCCConfig
)
from seqnado.config.user_input import load_genome_configs
from seqnado_web.design_utils import generate_design_from_fastqs

def scan_directory_for_fastqs(path: str) -> List[Dict[str, str]]:
    """
    Recursively scans a directory for FASTQ files.
    Returns a list of dicts with 'path' (absolute) and 'rel_path' (relative to search root).
    """
    root = Path(path).resolve()
    if not root.exists():
        return []

    fastqs = []
    # Patterns to look for
    patterns = ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"]
    
    # We use rglob for recursive search
    for pattern in patterns:
        for p in root.rglob(pattern):
            fastqs.append({
                "path": str(p),
                "rel_path": str(p.relative_to(root)),
                "name": p.name
            })
            
    # Sort by relative path
    fastqs.sort(key=lambda x: x["rel_path"])
    return fastqs

def create_project_structure(
    assay_name: str,
    project_name: str,
    genome_name: str,
    output_dir: str,
    selected_fastqs: List[str],
    config_options: Dict[str, str] = None
) -> str:
    """
    Creates the project directory structure, symlinks selected FASTQs,
    generates the design file from those symlinks, and creates the config.yaml.
    
    Args:
        assay_name: Name of the assay (rna, atac, etc.)
        project_name: Name of the project
        genome_name: Name of the genome config to use
        output_dir: Parent directory for the project
        selected_fastqs: List of absolute paths to FASTQ files to include
        config_options: Dict of extra config options from the wizard form
        
    Returns:
        Path to the generated project directory
    """
    if config_options is None:
        config_options = {}
    assay = Assay.from_clean_name(assay_name)
    
    # 1. Create Directories
    run_dir_name = f"{date.today().isoformat()}_{project_name}"
    run_dir = Path(output_dir) / run_dir_name
    fastq_dir = run_dir / "fastqs"
    fastq_dir.mkdir(parents=True, exist_ok=True)
    print(f"DEBUG: Created run_dir: {run_dir}")
    print(f"DEBUG: Created fastq_dir: {fastq_dir}")
    
    # 2. Symlink selected FASTQs
    count = 0
    for fq_path in selected_fastqs:
        src = Path(fq_path).resolve()
        if not src.exists():
            print(f"WARNING: Source file not found: {src}")
            continue
            
        # Create symlink in fastq_dir
        # We keep the filename
        dst = fastq_dir / src.name
        
        # Avoid duplicates or overwrite? 
        if dst.exists():
            print(f"DEBUG: Destination exists, skipping: {dst}")
            continue
            
        try:
            os.symlink(src, dst)
            count += 1
        except OSError as e:
            print(f"ERROR: Symlink failed for {src} -> {dst}: {e}")
            # Fallback to copy if symlink fails (e.g. cross-filesystem limitations)
            try:
                shutil.copy2(src, dst)
                count += 1
            except Exception as e2:
                 print(f"ERROR: Copy failed too: {e2}")

    print(f"DEBUG: Linked/Copied {count} files.")
            
    # 3. Generate Design from the *symlinked* files
    # We point to the new fastq_dir. The design generation logic handles pairing etc.
    design_filename = "metadata.csv"
    design_path_str = generate_design_from_fastqs(
        folder_path=str(fastq_dir),
        assay_name=assay_name,
        output_name=str((run_dir / design_filename).resolve())
    )
    design_path = Path(design_path_str)
    
    # 4. Create and Write Config
    # Create valid AssaySpecificConfig
    # For init, we used defaults which are generally None (disabled) or simple defaults.
    # The users can edit config.yaml later.
    
    # Load genome
    available_genomes = load_genome_configs(assay)
    genome_conf = available_genomes.get(genome_name)
    if not genome_conf:
         # Fallback: create a basic one with just name
         genome_conf = GenomeConfig(name=genome_name, index=None)

    project_conf = ProjectConfig(
        name=project_name,
        metadata_file=str(design_path.name)
    )

    # helper to get specific config class
    from seqnado.config.core import (
        ASSAY_CONFIG_MAP, ATACAssayConfig, ChIPAssayConfig, CATAssayConfig,
        RNAAssayConfig, MCCAssayConfig, SNPAssayConfig, MethylationAssayConfig, CRISPRAssayConfig
    )
    from seqnado import PileupMethod, PeakCallingMethod, QuantificationMethod
    
    # Build sub-configs from config_options
    bigwigs = None
    if config_options.get("make_bigwigs"):
        bigwigs = BigwigConfig(
            pileup_method=[PileupMethod(config_options.get("pileup_method", "deeptools"))],
            binsize=config_options.get("binsize", 10)
        )
    
    peak_calling = None
    if config_options.get("call_peaks"):
        peak_calling = PeakCallingConfig(
            method=[PeakCallingMethod(config_options.get("peak_calling_method", "lanceotron"))]
        )
    
    rna_quant = None
    if assay == Assay.RNA:
        rna_quant = RNAQuantificationConfig(
            method=QuantificationMethod(config_options.get("quant_method", "feature_counts")),
            run_deseq2=config_options.get("run_deseq2", False)
        )
    
    mcc_config = None
    if assay == Assay.MCC and config_options.get("mcc_viewpoints"):
        resolutions = [int(r.strip()) for r in config_options.get("mcc_resolutions", "100").split(",") if r.strip().isdigit()]
        mcc_config = MCCConfig(
            viewpoints=Path(config_options.get("mcc_viewpoints")),
            resolutions=resolutions if resolutions else [100]
        )
    
    # Build assay-specific config
    base_assay_cfg = {
        "bigwigs": bigwigs,
        "create_heatmaps": config_options.get("make_heatmaps", False),
    }
    
    if assay == Assay.ATAC:
        assay_config = ATACAssayConfig(
            **base_assay_cfg,
            tn5_shift=config_options.get("tn5_shift", True),
            peak_calling=peak_calling
        )
    elif assay == Assay.CHIP:
        assay_config = ChIPAssayConfig(
            **base_assay_cfg,
            peak_calling=peak_calling
        )
    elif assay == Assay.CAT:
        assay_config = CATAssayConfig(
            **base_assay_cfg,
            tn5_shift=config_options.get("tn5_shift", False),
            peak_calling=peak_calling
        )
    elif assay == Assay.RNA:
        assay_config = RNAAssayConfig(
            **base_assay_cfg,
            rna_quantification=rna_quant
        )
    elif assay == Assay.MCC:
        assay_config = MCCAssayConfig(
            **base_assay_cfg,
            mcc=mcc_config
        )
    else:
        ConfigClass = ASSAY_CONFIG_MAP.get(assay)
        assay_config = ConfigClass(**base_assay_cfg) if ConfigClass else None

    workflow_config = SeqnadoConfig(
        assay=assay,
        project=project_conf,
        genome=genome_conf,
        metadata=design_path.name,
        assay_config=assay_config
    )
    
    # Render config
    config_filename = f"config_{assay_name}.yaml"
    config_path = run_dir / config_filename
    
    try:
         from importlib.resources import files as _pkg_traversable
    except ImportError:
         from importlib.resources import files as _pkg_traversable

    tpl_trav = _pkg_traversable("seqnado.data").joinpath("config_template.jinja")
    
    with resources.as_file(tpl_trav) as tpl_path:
            render_config(
            template=Path(tpl_path),
            workflow_config=workflow_config,
            outfile=config_path,
            all_options=False
        )
        
    return str(run_dir)
