"""Project initialization and configuration helpers for pipeline tests."""

import glob
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

from .genome import fill_fastq_screen_config
from .utils import setup_genome_config


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
    genes_bed = test_data_path / "hg38_genes.bed"

    setup_genome_config(
        genome_config_file,
        star_index=resources.get("star_index"),
        bt2_index=resources["bt2_index"],
        chromsizes=resources["chromsizes"],
        gtf=resources["gtf"],
        blacklist=resources["blacklist"],
        genes_bed=genes_bed,
        fasta=resources.get("fasta")
        if any(x in assay.lower() for x in ["meth", "snp"])
        else None,
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
        assay = assay.replace("-rx", "")    

    # Generate config with proper flags
    result = subprocess.run(
        [
            "seqnado",
            "config",
            assay,
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
        str(run_directory / "**" / f"config_{assay}.yaml"), recursive=True
    )
    config_path = Path(config_files[0]) if config_files else None
    assert config_path.exists(), f"Config file not found for assay {assay}"

    if not config_path.exists():
        raise FileNotFoundError(
            f"Config file not found for assay {assay} at {config_path}"
        )

    with open(config_path) as f:
        config = yaml.safe_load(f)

    if config is None:
        raise ValueError(
            f"Loaded config for assay {assay} is None. Check the generated YAML file at {config_path}."
        )


    # Generate fastq_screen.conf for testing
    fastq_screen_config_path = run_directory / ".config" / "seqnado" / "fastq_screen.conf"
    fill_fastq_screen_config(fastq_screen_config_path, resources["bt2_index"])
    
    # Update config with test assay settings
    config["qc"]["run_fastq_screen"] = True
    config["genome"]["fastq_screen_config"] = str(fastq_screen_config_path)
    test_config_file = (
        Path(__file__).parent.parent / "assay_configs" / f"test_{assay}.yaml"
    )

    with open(test_config_file) as f:
        assay_config = yaml.safe_load(f)

    config.update(assay_config)

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



    result = subprocess.run(
        ["seqnado", "design", assay, "--no-interactive", "--accept-all-defaults"],
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
def multi_assay_run_directory(tmp_path_factory):
    """
    Fixture to provide a run directory for multi-assay tests.
    Returns a Path object to a unique temp directory.
    """
    return tmp_path_factory.mktemp("multi_assay_run")


def get_metadata_path(test_data_dir: Path, assay: str) -> Path:
    """Return the metadata CSV path for a given assay."""
    # Example: test_data_dir/metadata_atac.csv
    return test_data_dir / f"metadata_{assay}.csv"

