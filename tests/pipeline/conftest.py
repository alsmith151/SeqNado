from __future__ import annotations

import shutil
from pathlib import Path

import pytest

# All helpers are now in the helpers package
from helpers.config_utils import TestContext, make_test_paths
from helpers.data import ensure_fastqs_present
from helpers.genome import ensure_genome_resources
from helpers.project import (
    create_config_yaml,
    create_design_file,
    get_metadata_path,
    init_seqnado_project,
)
from helpers.utils import get_fastq_pattern


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
    make_test_paths(Path(__file__).resolve())

    return TestContext(pytestconfig, tmp_path_factory)


def ensure_seqnado_init(test_context, genome_resources, assay, monkeypatch):
    """
    Initialize seqnado project with required genome resources.
    """
    run_directory = test_context.run_directory(assay)
    resources = genome_resources(assay)  # Use genome_resources as a parameter

    # Ensure genome configuration is initialized
    init_seqnado_project(
        run_directory=run_directory,
        assay=assay,
        resources=resources,
        test_data_path=test_context.test_paths.test_data,
        monkeypatch=monkeypatch,
    )

    # Verify genome configuration file exists
    genome_config_file = run_directory / ".config" / "seqnado" / "genome_config.json"
    if not genome_config_file.exists():
        raise FileNotFoundError(f"Genome configuration file not found at {genome_config_file}")

    return resources


@pytest.fixture(scope="function")
def config_yaml_for_testing(
    test_context: TestContext,
    assay: str,
    monkeypatch: pytest.MonkeyPatch,
) -> Path:
    """Create and patch config YAML for testing."""
    run_directory = test_context.run_directory(assay)
    resources = genome_resources(assay)  # Use resources from genome_resources fixture
    config_path = create_config_yaml(run_directory, assay, monkeypatch, resources)
    return config_path


@pytest.fixture(scope="function")
def design(test_context: TestContext, assay: str, seqnado_run_dir: Path) -> Path:
    """Generate design file for the assay in the seqnado run directory."""
    assay_type = test_context.assay_type(assay)

    # Download FASTQ files
    ensure_fastqs_present(test_context.test_paths.fastq, [assay])
    fastq_source_dir = test_context.test_paths.fastq

    # Use the correct pattern from get_fastq_pattern
    pattern = get_fastq_pattern(assay_type)
    fastqs_to_copy = list(fastq_source_dir.glob(pattern))

    if not fastqs_to_copy:
        raise FileNotFoundError(
            f"No FASTQ files found for assay '{assay_type}' in {fastq_source_dir}"
        )

    # Move FASTQs to the run directory
    fastq_dest_dir = seqnado_run_dir / "fastqs"
    fastq_dest_dir.mkdir(parents=True, exist_ok=True)
    for fq in fastqs_to_copy:
        shutil.copy2(fq, fastq_dest_dir / fq.name)

    # Generate design file
    design_file = create_design_file(
        run_directory=seqnado_run_dir,
        assay=assay_type,
    )

    return design_file


@pytest.fixture(scope="function")
def multi_assay_configs(
    tmp_path_factory, test_context, monkeypatch, request, genome_resources
):
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

    # Initialize seqnado once in the shared run_dir for all assays
    # We need to collect all resources needed for all assays first
    all_resources = {}
    for assay in multi_assays:
        all_resources[assay] = genome_resources(assay)

    # Initialize with resources from the first assay (they should share common genome files)
    first_assay = multi_assays[0]
    init_seqnado_project(
        run_directory=run_dir,
        assay=first_assay,
        resources=all_resources[first_assay],
        test_data_path=test_context.test_paths.test_data,
        monkeypatch=monkeypatch,
    )

    for assay in multi_assays:
        # Get metadata and config paths
        metadata_path = get_metadata_path(test_context.test_paths.test_data, assay)

        # Create config YAML using helpers/project.py
        config_yaml = create_config_yaml(run_dir, assay, monkeypatch, all_resources[assay])

        configs[assay] = {"config": config_yaml, "metadata": metadata_path}
    return configs
