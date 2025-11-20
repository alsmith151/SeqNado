import os
import subprocess
from pathlib import Path

import pytest


@pytest.mark.pipeline
@pytest.mark.snakemake
@pytest.mark.requires_apptainer
@pytest.mark.slow
def test_pipeline(
    assay: str,
    assay_type: str,
    config_yaml_for_testing: Path,
    apptainer_args,
    cores: int,
    design: Path,
    test_profile_path: Path,
):
    """Execute the Snakemake workflow through the seqnado CLI using a generated config.

    Preconditions are handled by fixtures (init, config, design, test data, mounts).
    """
    res = subprocess.run(
        [
            "seqnado",
            "pipeline",
            assay_type,
            "-c",
            str(cores),
            "--configfile",
            str(config_yaml_for_testing),
            "--workflow-profile",
            str(test_profile_path),
        ],
        cwd=config_yaml_for_testing.parent,
        capture_output=True,
        text=True,
    )

    if res.returncode != 0:
        print("STDOUT:\n", res.stdout)
        print("STDERR:\n", res.stderr)
    assert res.returncode == 0

    # Check outputs relative to the test working directory
    test_dir = config_yaml_for_testing.parent
    assert not (test_dir / "seqnado_error.log").exists()
    assert (test_dir / "seqnado_output").exists()
    assert (test_dir / f"seqnado_output/{assay_type}/seqnado_report.html").exists()


@pytest.mark.pipeline
@pytest.mark.snakemake
@pytest.mark.requires_apptainer
@pytest.mark.slow
@pytest.mark.parametrize("multi_assays", [["atac", "chip", "rna"]])
def test_pipeline_multi(
    multi_assays: list[str],
    cores: int,
    test_profile_path: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """Execute the Snakefile_multi workflow with multiple assays."""
    # This test would need additional fixtures for multi-assay setup
    # For now, this is a placeholder that can be expanded
    pytest.skip("Multi-assay pipeline test not yet fully implemented")


def test_config_created(assay: str, config_yaml: Path, assay_type: str):
    assert os.path.exists(config_yaml), f"{assay_type} config file not created."


def test_design_created(assay: str, design: Path, assay_type: str):
    assert os.path.exists(design), f"{assay_type} design file not created."
