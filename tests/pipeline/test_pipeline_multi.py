"""Tests for multi-assay Snakefile_multi workflow.

These tests verify that multiple assays can be run together in a single
project using the Snakefile_multi orchestrator.
"""

import os
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
    cores: int,
    test_profile_path: Path,
    tmp_path: Path,
):
    """Execute the Snakefile_multi workflow with multiple assays.
    
    This test:
    1. Creates config and metadata files for each assay
    2. Runs Snakefile_multi to execute all assays
    3. Verifies outputs for each assay are created
    """
    pytest.skip("Multi-assay pipeline test requires additional fixtures - coming soon")
    
    # TODO: Implementation plan
    # 1. Create multi-assay fixture in conftest.py that:
    #    - Sets up configs for each assay in multi_assays
    #    - Creates metadata files for each assay
    #    - Returns a dict mapping assay -> config/metadata paths
    #
    # 2. Run snakemake with Snakefile_multi:
    #    res = subprocess.run([
    #        "snakemake",
    #        "-s", str(Path(__file__).parent.parent.parent / "seqnado/workflow/Snakefile_multi"),
    #        "-c", str(cores),
    #        "--workflow-profile", str(test_profile_path),
    #    ], cwd=test_dir, capture_output=True, text=True)
    #
    # 3. Verify outputs for each assay:
    #    for assay in multi_assays:
    #        assert (test_dir / f"seqnado_output/{assay}/seqnado_report.html").exists()
    #        assert (test_dir / f"seqnado_output/{assay}/logs/.complete").exists()
    #
    # 4. Verify summary file:
    #    assert (test_dir / "multi_assay_summary.txt").exists()
