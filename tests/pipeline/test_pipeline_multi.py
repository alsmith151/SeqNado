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
@pytest.mark.parametrize(
    "multi_assays",
    [
        [
            "atac",
            "chip",
            "meth",
            "rna",
            "snp",
        ]
    ],
)
def test_pipeline_multi(
    multi_assays: list[str],
    multi_assay_configs: dict[str, dict[str, Path]],
    multi_assay_run_directory: Path,
    cores: int,
    test_profile_path: Path,
):
    """Test running multiple assays in a single seqnado project using Snakefile_multi.
    Args:
        multi_assays: List of assay names to run together
        multi_assay_configs: Dict mapping assay names to their config and metadata paths
        multi_assay_run_directory: Path to the run directory for the multi-assay test
        cores: Number of cores to use for the pipeline
        test_profile_path: Path to the Snakemake profile for testing
    """
    # Find the Snakefile_multi
    snakefile_multi = (
        Path(__file__).parent.parent.parent / "seqnado" / "workflow" / "Snakefile_multi"
    )
    assert snakefile_multi.exists(), f"Snakefile_multi not found at {snakefile_multi}"

    res = subprocess.run(
        [
            "seqnado",
            "pipeline",
            "-c",
            str(cores),
            "--preset",
            "t",
        ],
        cwd=multi_assay_run_directory,
        capture_output=False,
        text=True,
    )

    assert res.returncode == 0, (
        f"Pipeline failed with return code {res.returncode}. See output above."
    )

    for assay in multi_assays:
        assert (multi_assay_run_directory / f"seqnado_output/{assay}").exists(), (
            f"Output directory not found for {assay}"
        )

        assert (
            multi_assay_run_directory / f"seqnado_output/{assay}/logs/.complete"
        ).exists(), f"No .complete log file found for {assay}"

    assert (multi_assay_run_directory / "multi_assay_summary.txt").exists(), (
        "multi_assay_summary.txt not found"
    )
