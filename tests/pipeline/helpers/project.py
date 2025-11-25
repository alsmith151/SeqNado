"""Project initialization and configuration helpers for pipeline tests."""
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
    genes_bed = test_data_path / "hg38_genes.bed"

    setup_genome_config(
        genome_config_file,
        star_index=resources.get("star_index"),
        bt2_index=resources["bt2_index"],
        chromsizes=resources["chromsizes"],
        gtf=resources["gtf"],
        blacklist=resources["blacklist"],
        genes_bed=genes_bed,
        fasta=resources.get("fasta") if any(x in assay.lower() for x in ["meth", "snp"]) else None,
    )

    assert genome_config_file.exists(), "genome_config.json not created"


def create_config_yaml(
    run_directory: Path,
    assay: str,
    monkeypatch: pytest.MonkeyPatch,
) -> Path:
    """
    Generate config YAML for the assay.

    Args:
        run_directory: Directory for the test run
        assay: Assay type
        monkeypatch: Pytest monkeypatch fixture

    Returns:
        Path to the generated config file
    """
    from datetime import datetime

    # Set environment for config generation
    monkeypatch.setenv("SEQNADO_CONFIG", str(run_directory))
    monkeypatch.setenv("HOME", str(run_directory))

    # Create project directory with datestamp (matching original behavior)
    date = datetime.now().strftime("%Y-%m-%d")
    project_dir = run_directory / f"{date}_{assay}_test"
    project_dir.mkdir(parents=True, exist_ok=True)

    # Generate config with proper flags
    config_path = project_dir / f"config_{assay}.yaml"
    result = subprocess.run(
        [
            "seqnado",
            "config",
            assay,
            "--no-interactive",
            "--no-make-dirs",
            "--render-options",
            "-o",
            str(config_path),
        ],
        cwd=run_directory,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"seqnado config failed:\nSTDERR: {result.stderr}\nSTDOUT: {result.stdout}"
    )
    assert config_path.exists(), f"Config file not created at {config_path}"

    return config_path


def patch_config_yaml(config_path: Path, assay: str) -> Path:
    """
    Patch the config YAML with test-specific settings.

    Args:
        config_path: Path to the config file
        assay: Assay type

    Returns:
        Path to the patched config file
    """
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Disable plotting for faster tests (matching original behavior)
    if "assay_config" in config:
        config["assay_config"]["plot_with_plotnado"] = False

    with open(config_path, "w") as f:
        yaml.dump(config, f, sort_keys=False)

    return config_path


def copy_fastqs(
    fastq_source_dir: Path,
    run_directory: Path,
    assay: str,
) -> list[Path]:
    """
    Copy FASTQ files for the assay directly to the run directory.

    Args:
        fastq_source_dir: Directory containing FASTQ files
        run_directory: Destination directory
        assay: Assay type

    Returns:
        List of copied FASTQ file paths
    """
    pattern = get_fastq_pattern(assay)
    fastqs_to_copy = list(fastq_source_dir.glob(pattern))

    if not fastqs_to_copy:
        raise FileNotFoundError(
            f"No FASTQ files matching pattern '{pattern}' for assay '{assay}' in {fastq_source_dir}"
        )

    # Copy FASTQs directly to run directory (not to a subdirectory)
    copied_fastqs = []
    for fq in fastqs_to_copy:
        dest = run_directory / fq.name
        if not dest.exists():
            shutil.copy2(fq, dest)
        copied_fastqs.append(dest)

    assert copied_fastqs, f"No FASTQ files copied for assay {assay} to {run_directory}"
    return sorted(copied_fastqs)


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
    design_file = run_directory / f"design_{assay}.csv"

    result = subprocess.run(
        ["seqnado", "design", assay, "--output", str(design_file)],
        cwd=run_directory,
        capture_output=True,
        text=True,
        input="y\nn\n",  # Accept first prompt, skip optional prompts
    )
    assert result.returncode == 0, (
        f"seqnado design failed:\nSTDERR: {result.stderr}\nSTDOUT: {result.stdout}"
    )
    assert design_file.exists(), f"Design file not created at {design_file}"

    return design_file
