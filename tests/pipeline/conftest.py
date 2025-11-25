from __future__ import annotations

from pathlib import Path

import pytest

# All helpers are now in the helpers package
from helpers.config_utils import TestContext, make_test_paths
from helpers.data import ensure_fastqs_present
from helpers.genome import ensure_genome_resources
from helpers.multi_assay import multi_assay_configs, multi_assay_run_directory
from helpers.project import (
    copy_fastqs,
    create_config_yaml,
    create_design_file,
    init_seqnado_project,
    patch_config_yaml,
)


# Fixture for all genome resources
@pytest.fixture(scope="session")
def genome_resources(test_context: TestContext) -> dict:
    """
    Provides a function to get genome resources for an assay.
    Returns a function that takes an assay and returns a dict of resource paths.
    """
    def _get_resources(assay: str) -> dict:
        return ensure_genome_resources(test_context.test_paths.genome, assay)
    return _get_resources


@pytest.fixture(scope="function")
def seqnado_run_dir(config_yaml_for_testing: Path) -> Path:
    """Return the seqnado run directory."""
    return Path(config_yaml_for_testing).parent


@pytest.fixture(scope="session")
def test_profile_path(test_context: TestContext) -> Path:
    """Return the path to the test profile for Snakemake."""
    return test_context.test_paths.test_profile


@pytest.fixture(scope="session")
def cores(test_context: TestContext) -> int:
    """Return the number of cores to use for pipeline tests."""
    return test_context.cores


@pytest.fixture(scope="session")
def test_context(pytestconfig, tmp_path_factory) -> TestContext:
    """Session-scoped test context with paths and settings."""

    selected_assays = pytestconfig.getoption("--assays")
    if isinstance(selected_assays, str):
        selected_assays = [a.strip() for a in selected_assays.split(",") if a.strip()]
    elif not selected_assays:
        selected_assays = []
    test_path = make_test_paths(Path(__file__).resolve())
    # Ensure FASTQ files are present for selected assays
    ensure_fastqs_present(test_path.test_data, selected_assays)
    return TestContext(pytestconfig, tmp_path_factory)


@pytest.fixture(scope="function")
def ensure_seqnado_init(
    test_context: TestContext,
    genome_resources,
    assay: str,
    monkeypatch: pytest.MonkeyPatch,
):
    """
    Initialize seqnado project with required genome resources.
    """
    run_directory = test_context.run_directory(assay)
    resources = genome_resources(assay)

    init_seqnado_project(
        run_directory=run_directory,
        assay=assay,
        resources=resources,
        test_data_path=test_context.test_paths.test_data,
        monkeypatch=monkeypatch,
    )


@pytest.fixture(scope="function")
def config_yaml_for_testing(
    test_context: TestContext,
    assay: str,
    monkeypatch: pytest.MonkeyPatch,
    ensure_seqnado_init,
) -> Path:
    """Create and patch config YAML for testing."""
    run_directory = test_context.run_directory(assay)
    config_path = create_config_yaml(run_directory, assay, monkeypatch)
    return patch_config_yaml(config_path, assay)


@pytest.fixture(scope="function")
def design(test_context: TestContext, assay: str, seqnado_run_dir: Path) -> Path:
    """Generate design file for the assay in the seqnado run directory."""
    assay_type = test_context.assay_type(assay)

    # Copy FASTQ files to run directory
    copy_fastqs(
        fastq_source_dir=test_context.test_paths.fastq,
        run_directory=seqnado_run_dir,
        assay=assay_type,
    )

    # Generate design file
    design_file = create_design_file(
        run_directory=seqnado_run_dir,
        assay=assay_type,
    )

    return design_file
