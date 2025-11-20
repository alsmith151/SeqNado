"""Tests for multi-assay Snakefile_multi workflow.

These tests verify that multiple assays can be run together in a single
project using the Snakefile_multi orchestrator.
"""

import subprocess
from pathlib import Path

import pytest


@pytest.mark.pipeline
@pytest.mark.snakemake
@pytest.mark.requires_apptainer
@pytest.mark.slow
@pytest.mark.parametrize("multi_assays", [["atac", "chip"]])
def test_pipeline_multi(
    multi_assays: list[str],
    multi_assay_configs: dict[str, dict[str, Path]],
    multi_assay_run_directory: Path,
    cores: int,
    test_profile_path: Path,
):
    """Execute the Snakefile_multi workflow with multiple assays.
    
    This test:
    1. Uses multi_assay_configs fixture to set up configs/metadata for each assay
    2. Runs Snakefile_multi to execute all assays
    3. Verifies outputs for each assay are created
    """
    # Find the Snakefile_multi
    snakefile_multi = Path(__file__).parent.parent.parent / "seqnado" / "workflow" / "Snakefile_multi"
    assert snakefile_multi.exists(), f"Snakefile_multi not found at {snakefile_multi}"
    
    # Run the multi-assay pipeline
    res = subprocess.run(
        [
            "snakemake",
            "-s", str(snakefile_multi),
            "-c", str(cores),
            "--workflow-profile", str(test_profile_path),
        ],
        cwd=multi_assay_run_directory,
        capture_output=True,
        text=True,
    )
    
    if res.returncode != 0:
        print("STDOUT:\n", res.stdout)
        print("STDERR:\n", res.stderr)
    assert res.returncode == 0, "Snakefile_multi failed"
    
    # Verify outputs for each assay
    for assay in multi_assays:
        assert (multi_assay_run_directory / f"seqnado_output/{assay}/seqnado_report.html").exists(), \
            f"Report not found for {assay}"
        assert (multi_assay_run_directory / f"seqnado_output/{assay}/logs/.complete").exists(), \
            f"Completion marker not found for {assay}"
    
    # Verify summary file is created
    assert (multi_assay_run_directory / "multi_assay_summary.txt").exists(), \
        "Multi-assay summary not created"

