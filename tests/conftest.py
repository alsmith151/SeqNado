"""Shared pytest configuration and fixtures for SeqNado tests.

This file defines global pytest options, markers and lightweight fixtures used
across unit, CLI and pipeline tests. Heavyweight pipeline fixtures live in
``tests/pipeline/conftest.py`` to keep scope clear.
"""
from pathlib import Path
import shutil
import pytest


@pytest.fixture(scope="session")
def test_data_dir():
    """Provide path to test data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def test_fastq_dir(test_data_dir):
    """Provide path to test FASTQ directory."""
    return test_data_dir / "fastq"


@pytest.fixture
def sample_metadata_dict():
    """Provide sample metadata dictionary for testing."""
    return {
        "scaling_group": "treatment",
        "condition": "high_dose",
        "spikein_ratio": 0.25
    }


@pytest.fixture
def multiple_metadata_dicts():
    """Provide multiple metadata dictionaries for testing."""
    return [
        {"scaling_group": "control", "condition": "baseline"},
        {"scaling_group": "treatment", "condition": "low_dose", "spikein_ratio": 0.1},
        {"scaling_group": "treatment", "condition": "high_dose", "spikein_ratio": 0.25},
    ]


@pytest.fixture
def sample_names():
    """Provide sample names for testing."""
    return ["sample1", "sample2", "sample3", "control1", "control2"]


def pytest_configure(config):
    """Configure pytest with custom markers and settings."""
    config.addinivalue_line(
        "markers", "unit: mark test as a unit test"
    )
    config.addinivalue_line(
        "markers", "integration: mark test as an integration test"
    )
    config.addinivalue_line(
        "markers", "edge_case: mark test as an edge case test"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )
    config.addinivalue_line(
        "markers", "requires_data: mark test as requiring data files"
    )
    config.addinivalue_line(
        "markers", "pipeline: end-to-end Snakemake pipeline tests (opt-in)"
    )
    config.addinivalue_line(
        "markers", "snakemake: tests that invoke Snakemake via subprocess"
    )
    config.addinivalue_line(
        "markers", "requires_apptainer: tests that require Apptainer/Singularity"
    )


def pytest_addoption(parser: pytest.Parser) -> None:
    """CLI options to control pipeline tests.

    --run-pipeline: Opt-in flag to execute slow pipeline tests.
    --assays: Comma-separated list of assays to run for pipeline tests.
    --cores: Number of cores to use when running the pipeline.
    """
    group = parser.getgroup("seqnado")
    group.addoption(
        "--run-pipeline",
        action="store_true",
        default=False,
        help="Run slow pipeline tests (disabled by default)",
    )
    group.addoption(
        "--assays",
        action="store",
        default="chip",
        help=(
            "Comma-separated list of assays to test (choices: atac, chip, chip-rx, "
            "rna, rna-rx, snp, cat, meth, mcc). Default: chip"
        ),
    )
    group.addoption(
        "--cores",
        action="store",
        type=int,
        default=2,
        help="Number of CPU cores to allocate to pipeline tests (default: 2)",
    )


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]):
    """Modify test collection to add markers based on test names and locations."""
    # Determine environment/tool availability once
    apptainer_available = bool(shutil.which("apptainer") or shutil.which("singularity"))

    for item in items:
        # Add unit marker to all tests by default
        if not any(mark.name in ["integration", "slow"] for mark in item.iter_markers()):
            item.add_marker(pytest.mark.unit)
        
        # Add edge_case marker to tests with "edge" in name
        if "edge" in item.name.lower():
            item.add_marker(pytest.mark.edge_case)
        
        # Add requires_data marker to tests that use test data
        if any(fixture in item.fixturenames for fixture in ["test_data_dir", "test_fastq_dir"]):
            item.add_marker(pytest.mark.requires_data)

        # If tests are marked as pipeline and user did not opt in, skip them
        if any(m.name == "pipeline" for m in item.iter_markers()):
            if not config.getoption("--run-pipeline"):
                item.add_marker(pytest.mark.skip(reason="Use --run-pipeline to enable pipeline tests"))

        # If a test requires Apptainer and it's not available, skip it
        if any(m.name == "requires_apptainer" for m in item.iter_markers()):
            if not apptainer_available:
                item.add_marker(pytest.mark.skip(reason="Apptainer/Singularity not found in PATH"))


# Custom assertions for SeqNado
def assert_valid_sample_name(sample_name: str):
    """Assert that a sample name is valid."""
    assert isinstance(sample_name, str), "Sample name must be a string"
    assert len(sample_name) > 0, "Sample name cannot be empty"
    assert not sample_name.isspace(), "Sample name cannot be only whitespace"


def assert_valid_file_path(path: Path):
    """Assert that a file path is valid and exists."""
    assert isinstance(path, Path), "Path must be a Path object"
    assert path.exists(), f"Path does not exist: {path}"
    assert path.is_file(), f"Path is not a file: {path}"


# Make custom assertions available to all tests
pytest.assert_valid_sample_name = assert_valid_sample_name
pytest.assert_valid_file_path = assert_valid_file_path


# -------------------------
# Global lightweight fixtures
# -------------------------

@pytest.fixture(scope="session")
def cores(pytestconfig: pytest.Config) -> int:
    """Number of cores to use for pipeline executions."""
    return int(pytestconfig.getoption("--cores"))


@pytest.fixture(scope="session")
def selected_assays(pytestconfig: pytest.Config) -> list[str]:
    """List of assays to parametrize pipeline tests with."""
    assays_str: str = pytestconfig.getoption("--assays") or "chip"
    assays = [a.strip() for a in assays_str.split(",") if a.strip()]
    # Normalize to known set; leave as-is to allow custom experimental names if needed
    return assays or ["chip"]


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """Parametrize tests that accept the 'assay' fixture using --assays values."""
    if "assay" in metafunc.fixturenames:
        config = metafunc.config
        assays_str: str = config.getoption("--assays") or "chip"
        assays = [a.strip() for a in assays_str.split(",") if a.strip()]
        metafunc.parametrize("assay", assays or ["chip"])
