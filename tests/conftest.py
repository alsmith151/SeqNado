"""Shared pytest configuration and fixtures for SeqNado tests."""
import pathlib
import pytest


@pytest.fixture(scope="session")
def test_data_dir():
    """Provide path to test data directory."""
    return pathlib.Path(__file__).parent / "data"


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


def pytest_collection_modifyitems(config, items):
    """Modify test collection to add markers based on test names and locations."""
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


# Custom assertions for SeqNado
def assert_valid_sample_name(sample_name: str):
    """Assert that a sample name is valid."""
    assert isinstance(sample_name, str), "Sample name must be a string"
    assert len(sample_name) > 0, "Sample name cannot be empty"
    assert not sample_name.isspace(), "Sample name cannot be only whitespace"


def assert_valid_file_path(path: pathlib.Path):
    """Assert that a file path is valid and exists."""
    assert isinstance(path, pathlib.Path), "Path must be a pathlib.Path object"
    assert path.exists(), f"Path does not exist: {path}"
    assert path.is_file(), f"Path is not a file: {path}"


# Make custom assertions available to all tests
pytest.assert_valid_sample_name = assert_valid_sample_name
pytest.assert_valid_file_path = assert_valid_file_path
