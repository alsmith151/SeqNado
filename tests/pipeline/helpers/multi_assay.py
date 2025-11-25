from pathlib import Path

import pytest

from .data import ensure_fastqs_present
from .project import create_config_yaml


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


@pytest.fixture(scope="function")
def multi_assay_configs(tmp_path_factory, test_context, monkeypatch, request):
    """
    For a list of assays, set up config and metadata files for each assay.
    Returns: {assay: {"config": config_path, "metadata": metadata_path}}
    """
    # Get the list of assays from the test parameterization or request
    multi_assays = getattr(request, "param", None)
    if multi_assays is None and hasattr(request, "node"):
        # Try to get from test function arguments
        multi_assays = request.node.funcargs.get("multi_assays", None)
    if multi_assays is None:
        # Fallback: use a default
        multi_assays = ["atac", "rna"]
    # Ensure FASTQ files are present for all requested assays
    ensure_fastqs_present(test_context.test_paths.fastq, multi_assays)
    # Set up a run directory for the multi-assay test
    run_dir = tmp_path_factory.mktemp("multi_assay_run")
    configs = {}
    for assay in multi_assays:
        # Get metadata and config paths
        metadata_path = get_metadata_path(test_context.test_paths.test_data, assay)

        # Create config YAML using helpers/project.py
        config_yaml = create_config_yaml(run_dir, assay, monkeypatch)

        configs[assay] = {"config": config_yaml, "metadata": metadata_path}
    return configs
