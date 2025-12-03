import os
import subprocess
from pathlib import Path

import pytest


@pytest.mark.pipeline
@pytest.mark.snakemake
@pytest.mark.requires_apptainer
@pytest.mark.slow
def test_pipeline(
    assay,
    test_context,
    config_yaml_for_testing: Path,
    design: Path,
    test_profile_path: Path,
):
    assay_type = test_context.assay_type(assay)
    cores = test_context.cores
    res = subprocess.run(
        [
            "seqnado",
            "pipeline",
            assay_type,
            "-c",
            str(cores),
            "--configfile",
            str(config_yaml_for_testing),
            "--preset",
            "t",
        ],
        cwd=config_yaml_for_testing.parent,
        capture_output=True,
        text=True,
    )

    # Print output for debugging if the test fails
    if res.returncode != 0:
        print("\n=== STDOUT ===")
        print(res.stdout)
        print("\n=== STDERR ===")
        print(res.stderr)

    assert res.returncode == 0, f"Pipeline failed with return code {res.returncode}. See output above."
    test_dir = config_yaml_for_testing.parent
    assert not (test_dir / "seqnado_error.log").exists()
    assert (test_dir / "seqnado_output").exists()
    assert (test_dir / f"seqnado_output/{assay_type}/seqnado_report.html").exists()
