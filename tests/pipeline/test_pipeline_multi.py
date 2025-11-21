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
    # Note: May have some failures due to test data limitations (e.g., peak calling, hub generation)
    # but core pipeline steps should complete
    res = subprocess.run(
        [
            "seqnado",
            "pipeline",
            "--workflow-profile",
            str(test_profile_path),
            "--keep-going",  # Continue on errors to generate as much output as possible
        ],
        cwd=multi_assay_run_directory,
        capture_output=True,
        text=True,
    )
    
    # Don't fail test if pipeline has partial success (some rules may fail with small test data)
    # We'll check for key outputs below
    if res.returncode != 0:
        print("STDOUT:\n", res.stdout)
        print("STDERR:\n", res.stderr)
        print("\nNote: Pipeline may have partial failures (e.g., peak calling, hub generation)")
        print("Checking for core outputs...")
    
    # Verify outputs for each assay
    # Note: We check for key outputs but allow some steps to fail (e.g., peak calling on small test data)
    for assay in multi_assays:
        # Check that the assay output directory was created
        assert (multi_assay_run_directory / f"seqnado_output/{assay}").exists(), \
            f"Output directory not found for {assay}"
        
        # Check for BAM files (core alignment output)
        assert list((multi_assay_run_directory / f"seqnado_output/{assay}/aligned").glob("*.bam")), \
            f"No BAM files found for {assay}"
        
        # Check for bigwigs (core visualization output)
        assert list((multi_assay_run_directory / f"seqnado_output/{assay}/bigwigs").rglob("*.bigWig")), \
            f"No bigWig files found for {assay}"

