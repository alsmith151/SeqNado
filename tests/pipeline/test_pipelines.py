import os
import subprocess
from pathlib import Path

import pytest


@pytest.mark.pipeline
def test_config(assay, test_context, config_yaml_for_testing: Path):
    assay_type = test_context.assay_type(assay)
    assert os.path.exists(config_yaml_for_testing), (
        f"{assay_type} config file not created."
    )


@pytest.mark.pipeline
def test_design(assay, test_context, design: Path):
    assay_type = test_context.assay_type(assay)
    assert os.path.exists(design), f"{assay_type} design file not created."


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
            "ls",
        ],
        cwd=config_yaml_for_testing.parent,
        capture_output=False,
        text=True,
    )
    assert res.returncode == 0
    test_dir = config_yaml_for_testing.parent
    assert not (test_dir / "seqnado_error.log").exists()
    assert (test_dir / "seqnado_output").exists()
    assert (test_dir / f"seqnado_output/{assay_type}/seqnado_report.html").exists()
