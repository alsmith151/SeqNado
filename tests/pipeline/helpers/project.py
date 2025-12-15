import glob
import json
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

from .utils import get_fastq_pattern, setup_genome_config


def init_seqnado_project(
    run_directory: Path,
    assay: str,
    resources: dict,
    test_data_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """
    Initialize a seqnado project for testing.

    Args:
        run_directory: Directory for the test run
        assay: Assay type (e.g., 'chip', 'rna', 'mcc')
        resources: Dict of genome resources from ensure_genome_resources()
        test_data_path: Base path for test data
        monkeypatch: Pytest monkeypatch fixture for environment variables

    Sets up environment, copies necessary files, runs seqnado init, and writes genome config.
    """

    # Set environment variables
    monkeypatch.setenv("SEQNADO_CONFIG", str(run_directory))
    monkeypatch.setenv("HOME", str(run_directory))

    # MCC-specific setup
    if assay.lower() == "mcc":
        viewpoints_path = test_data_path / "genome" / "mcc_viewpoints.bed"
        monkeypatch.setenv("SEQNADO_MCC_VIEWPOINTS", str(viewpoints_path))

        # Copy viewpoints file if available
        if "viewpoints" in resources:
            dest = run_directory / "mcc_viewpoints.bed"
            src = resources["viewpoints"]
            if not dest.exists() and src.exists():
                shutil.copy2(src, dest)

    # Run seqnado init
    result = subprocess.run(
        ["seqnado", "init"],
        input="y\n",
        text=True,
        cwd=run_directory,
        capture_output=True,
    )
    assert result.returncode == 0, f"seqnado init failed: {result.stderr}"

    # Write genome config
    genome_config_file = run_directory / ".config" / "seqnado" / "genome_config.json"

    # Copy hg38_genes.bed from tests/data to test_output/data if it doesn't exist
    genes_bed_source = test_data_path.parent.parent / "tests" / "data" / "hg38_genes.bed"
    genes_bed = test_data_path / "hg38_genes.bed"
    if not genes_bed.exists() and genes_bed_source.exists():
        shutil.copy2(genes_bed_source, genes_bed)

    # Copy plotting_coordinates.bed from tests/data to test_output/data if it doesn't exist
    plot_coords_source = test_data_path.parent.parent / "tests" / "data" / "plotting_coordinates.bed"
    plot_coords = test_data_path / "plotting_coordinates.bed"
    if not plot_coords.exists() and plot_coords_source.exists():
        shutil.copy2(plot_coords_source, plot_coords)

    setup_genome_config(
        genome_config_file,
        star_index=resources.get("star_index"),
        bt2_index=resources["bt2_index"],
        chromsizes=resources["chromsizes"],
        gtf=resources["gtf"],
        blacklist=resources["blacklist"],
        genes_bed=genes_bed,
        fasta=resources.get("fasta"),
        assay=assay,
    )

    assert genome_config_file.exists(), "genome_config.json not created"


def create_config_yaml(
    run_directory: Path,
    assay: str,
    monkeypatch: pytest.MonkeyPatch,
    resources: dict,  # Added resources parameter
) -> Path:
    """
    Generate config YAML for the assay.

    Args:
        run_directory: Directory for the test run
        assay: Assay type
        monkeypatch: Pytest monkeypatch fixture
        resources: Dictionary containing genome resources

    Returns:
        Path to the generated config file
    """

    # Set environment for config generation
    monkeypatch.setenv("SEQNADO_CONFIG", str(run_directory))
    monkeypatch.setenv("HOME", str(run_directory))

    # amend -rx assays to base assay names
    if assay.endswith("-rx"):
        seq_assay = assay.replace("-rx", "")
    else:
        seq_assay = assay

    # Generate config with proper flags
    result = subprocess.run(
        [
            "seqnado",
            "config",
            seq_assay,
            "--no-interactive",
        ],
        cwd=run_directory,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"seqnado config failed:\nSTDERR: {result.stderr}\nSTDOUT: {result.stdout}"
    )

    # Find the generated config file (search recursively)
    config_files = glob.glob(
        str(run_directory / "**" / f"config_{seq_assay}.yaml"), recursive=True
    )
    config_path = Path(config_files[0]) if config_files else None
    assert config_path.exists(), f"Config file not found for assay {seq_assay}"

    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Load existing genome_config.json and get the genome entry for this assay
    genome_config_file = run_directory / ".config" / "seqnado" / "genome_config.json"
    with open(genome_config_file) as f:
        genome_config_data = json.load(f)

    # Determine which genome config to use
    if assay in genome_config_data:
        genome_key = assay
        genome_config = genome_config_data[assay]
    else:
        genome_key = "hg38"
        genome_config = genome_config_data.get("hg38", {})

    # Update genome entries in config with paths from genome_config.json
    config["genome"]["name"] = genome_key

    # Update the index prefix based on assay type
    if "rna" in assay.lower() and genome_config.get("star_index"):
        config["genome"]["index"]["prefix"] = genome_config.get("star_index")
        config["genome"]["index"]["type"] = "STAR"
    elif genome_config.get("bt2_index"):
        config["genome"]["index"]["prefix"] = genome_config.get("bt2_index")
        config["genome"]["index"]["type"] = "Bowtie2"
    else:
        config["genome"]["index"]["prefix"] = None
        config["genome"]["index"]["type"] = None
        

    # Update other genome fields (chromosome_sizes from JSON maps to both field names)
    config["genome"]["chromosome_sizes"] = genome_config.get("chromosome_sizes")
    config["genome"]["gtf"] = genome_config.get("gtf")
    config["genome"]["blacklist"] = genome_config.get("blacklist")
    config["genome"]["organism"] = "Homo sapiens"
    if "fasta" in genome_config:
        config["genome"]["fasta"] = genome_config.get("fasta")

    # Generate fastq_screen.conf for testing
    genome_path = Path(resources["bt2_index"]).parent.parent
    fastq_screen_config_path = genome_path / "fastq_screen.conf"

    # Update config with test assay settings
    config["qc"]["run_fastq_screen"] = True
    config["genome"]["fastq_screen_config"] = str(fastq_screen_config_path)

    # Also update the third_party_tools config which is what the Snakemake rule uses
    if "third_party_tools" in config and "fastq_screen" in config["third_party_tools"]:
        config["third_party_tools"]["fastq_screen"]["config"] = str(
            fastq_screen_config_path
        )
    test_config_file = (
        Path(__file__).parent.parent / "assay_configs" / f"test_{assay}.yaml"
    )

    with open(test_config_file) as f:
        assay_config = yaml.safe_load(f)

    config.update(assay_config)

    # Add genome directory bind mount for Singularity/Apptainer
    # This allows the container to access genome files (e.g., bt2 indexes for fastq_screen)
    # We need to update the Snakemake profile config, not the workflow config
    genome_dir = Path(resources["bt2_index"]).parent.parent.resolve()

    # Also mount the parent test_data directory so plotting_coordinates.bed is accessible
    test_data_dir = genome_dir.parent.resolve()

    # Find and update the test profile configuration
    import site
    import sys

    # Try to find the seqnado package location
    seqnado_paths = [p for p in sys.path if 'seqnado' in p and 'site-packages' in p]
    if not seqnado_paths:
        # Fall back to searching site-packages
        for site_pkg in site.getsitepackages():
            test_profile_config = Path(site_pkg) / "seqnado" / "workflow" / "envs" / "profiles" / "profile_test" / "config.v8+.yaml"
            if test_profile_config.exists():
                break
    else:
        test_profile_config = Path(seqnado_paths[0]) / "seqnado" / "workflow" / "envs" / "profiles" / "profile_test" / "config.v8+.yaml"

    # Update the profile config with apptainer-args
    if test_profile_config.exists():
        with open(test_profile_config) as f:
            profile_config = yaml.safe_load(f)

        # Bind mount test_data_dir instead of just genome_dir to ensure plotting_coordinates.bed is accessible
        bind_arg = f"--bind {test_data_dir}:{test_data_dir}"
        if "apptainer-args" in profile_config:
            # Check if this bind mount is already present
            if str(test_data_dir) not in profile_config["apptainer-args"]:
                profile_config["apptainer-args"] += f" {bind_arg}"
        else:
            profile_config["apptainer-args"] = bind_arg

        with open(test_profile_config, "w") as f:
            yaml.dump(profile_config, f, sort_keys=False)

    # Fix plotting coordinates path to use test_output/data instead of package directory
    if "assay_config" in config and "plotting" in config["assay_config"]:
        # Get the test data directory from the genome config path
        # chromosome_sizes is in test_output/data/genome/, so go up one level to get test_output/data/
        test_data_dir = Path(genome_config.get("chromosome_sizes")).parent.parent
        plot_coords = test_data_dir / "plotting_coordinates.bed"
        # Always update the path to point to test_output/data, regardless of whether it exists yet
        # The conftest fixture will ensure the file is copied before the test runs
        config["assay_config"]["plotting"]["coordinates"] = str(plot_coords)

    with open(config_path, "w") as f:
        yaml.dump(config, f, sort_keys=False)

    return config_path


def create_design_file(
    run_directory: Path,
    assay: str,
) -> Path:
    """
    Generate design CSV file for the assay.

    Args:
        run_directory: Directory for the test run
        assay: Assay type

    Returns:
        Path to the generated design file
    """

    # amend -rx assays to base assay names
    if assay.endswith("-rx"):
        assay = assay.replace("-rx", "")

    # Expand the glob pattern to get actual FASTQ file paths
    pattern = get_fastq_pattern(assay)

    # For multiomics, FASTQs are in fastqs/{assay}/ subdirectory
    fastq_dir = run_directory / "fastqs" / assay

    # Fall back to fastqs/ if assay-specific directory doesn't exist
    if not fastq_dir.exists():
        fastq_dir = run_directory / "fastqs"

    fastq_files = sorted(fastq_dir.glob(pattern))

    if not fastq_files:
        raise FileNotFoundError(
            f"No FASTQ files found matching pattern '{pattern}' in {fastq_dir}"
        )

    # Build command with expanded file paths
    cmd = [
        "seqnado",
        "design",
        assay,
        "--no-interactive",
        "--accept-all-defaults",
    ]

    # Add relative paths to FASTQ files
    # For multiomics structure: fastqs/{assay}/{file}.fastq.gz
    if (run_directory / "fastqs" / assay).exists():
        cmd.extend([f"fastqs/{assay}/{f.name}" for f in fastq_files])
    else:
        # Fall back to old structure: fastqs/{file}.fastq.gz
        cmd.extend([f"fastqs/{f.name}" for f in fastq_files])

    result = subprocess.run(
        cmd,
        cwd=run_directory,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"seqnado design failed:\nSTDERR: {result.stderr}\nSTDOUT: {result.stdout}"
    )

    design_file = run_directory / f"metadata_{assay}.csv"
    assert design_file.exists(), f"Design file not created at {design_file}"

    return design_file


@pytest.fixture(scope="function")
def multiomics_run_directory(tmp_path_factory):
    """
    Fixture to provide a run directory for multi-assay tests.
    Returns a Path object to a unique temp directory.
    """
    return tmp_path_factory.mktemp("multiomics_run")


def get_metadata_path(test_data_dir: Path, assay: str) -> Path:
    """Return the metadata CSV path for a given assay."""
    # Example: test_data_dir/metadata_atac.csv
    return test_data_dir / f"metadata_{assay}.csv"
