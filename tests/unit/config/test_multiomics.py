"""Unit tests for seqnado.config.multiomics module."""

from pathlib import Path

import pytest

from seqnado import Assay
from seqnado.outputs.multiomics import (
    find_assay_config_paths,
    find_metadata_paths,
    validate_config_and_metadata,
    none_str_to_none,
)

from seqnado.config import MultiomicsConfig

class TestNoneStrToNone:
    """Tests for the none_str_to_none validator function."""

    def test_converts_string_none_to_none(self):
        """Test that string 'none' is converted to None."""
        assert none_str_to_none("none") is None

    def test_converts_uppercase_none_to_none(self):
        """Test that uppercase 'NONE' is converted to None."""
        assert none_str_to_none("NONE") is None

    def test_converts_mixed_case_none_to_none(self):
        """Test that mixed case 'NoNe' is converted to None."""
        assert none_str_to_none("NoNe") is None

    def test_converts_none_with_whitespace_to_none(self):
        """Test that 'none' with whitespace is converted to None."""
        assert none_str_to_none("  none  ") is None
        assert none_str_to_none("\tnone\n") is None

    def test_returns_non_none_strings_unchanged(self):
        """Test that non-'none' strings are returned unchanged."""
        assert none_str_to_none("something") == "something"
        assert none_str_to_none("") == ""
        assert none_str_to_none("not_none") == "not_none"

    def test_returns_non_string_values_unchanged(self):
        """Test that non-string values are returned unchanged."""
        assert none_str_to_none(None) is None
        assert none_str_to_none(42) == 42
        assert none_str_to_none([1, 2, 3]) == [1, 2, 3]
        assert none_str_to_none({"key": "value"}) == {"key": "value"}


class TestFindAssayConfigs:
    """Tests for the find_assay_configs function."""

    def test_finds_single_assay_config(self, tmp_path: Path):
        """Test finding a single assay config and metadata pair."""
        # Create config and metadata files
        config_file = tmp_path / "config_chip.yaml"
        config_file.write_text("test: config")

        metadata_file = tmp_path / "metadata_chip.csv"
        metadata_file.write_text("sample,condition\n")

        config_files = find_assay_config_paths(tmp_path)
        metadata_files = find_metadata_paths(tmp_path)

        assert len(config_files) == 1
        assert Assay.CHIP in config_files
        assert config_files[Assay.CHIP] == Path(config_file)
        
        assert len(metadata_files) == 1
        assert Assay.CHIP in metadata_files
        assert metadata_files[Assay.CHIP] == Path(metadata_file)

    def test_finds_multiple_assay_configs(self, tmp_path: Path):
        """Test finding multiple assay configs and metadata pairs."""
        assays = ["chip", "rna", "atac"]

        for assay in assays:
            config_file = tmp_path / f"config_{assay}.yaml"
            config_file.write_text(f"assay: {assay}")

            metadata_file = tmp_path / f"metadata_{assay}.csv"
            metadata_file.write_text("sample,condition\n")

        config_files = find_assay_config_paths(tmp_path)
        metadata_files = find_metadata_paths(tmp_path)

        assert len(config_files) == 3
        assert len(metadata_files) == 3

        for assay_name in assays:
            assay = Assay.from_clean_name(assay_name)
            assert assay in config_files
            assert assay in metadata_files
            assert config_files[assay] == tmp_path / f"config_{assay_name}.yaml"
            assert metadata_files[assay] == tmp_path / f"metadata_{assay_name}.csv"

    def test_raises_error_when_metadata_missing(self, tmp_path: Path):
        """Test that FileNotFoundError is raised when metadata file is missing."""
        # Create config file but not metadata
        config_file = tmp_path / "config_chip.yaml"
        config_file.write_text("test: config")

        config_files = find_assay_config_paths(tmp_path)
        metadata_files = find_metadata_paths(tmp_path)

        with pytest.raises(FileNotFoundError) as exc_info:
            validate_config_and_metadata(config_files, metadata_files)

        assert "Missing metadata file" in str(exc_info.value)
        assert "CHIP" in str(exc_info.value)

    def test_empty_directory_returns_empty_dicts(self, tmp_path: Path):
        """Test that empty directory returns empty dictionaries."""
        config_files = find_assay_config_paths(tmp_path)
        metadata_files = find_metadata_paths(tmp_path)

        assert config_files == {}
        assert metadata_files == {}

    def test_ignores_non_config_yaml_files(self, tmp_path: Path):
        """Test that non-config YAML files are ignored."""
        # Create files that don't match the config_*.yaml pattern
        (tmp_path / "other.yaml").write_text("test: data")
        (tmp_path / "settings.yaml").write_text("test: data")
        (tmp_path / "config.yaml").write_text("test: data")  # Missing underscore

        config_files = find_assay_config_paths(tmp_path)
        metadata_files = find_metadata_paths(tmp_path)

        assert config_files == {}
        assert metadata_files == {}

    def test_raises_error_for_unknown_assay_names(self, tmp_path: Path):
        """Test that unknown assay names raise ValueError."""
        config_file = tmp_path / "config_chip-rx.yaml"
        config_file.write_text("test: config")

        with pytest.raises(ValueError) as exc_info:
            find_assay_config_paths(tmp_path)
        
        assert "Unknown clean name: chip-rx" in str(exc_info.value)

    def test_partial_configs_raise_error(self, tmp_path: Path):
        """Test that having some configs with metadata and some without raises error."""
        # Create first pair correctly
        (tmp_path / "config_chip.yaml").write_text("test: config")
        (tmp_path / "metadata_chip.csv").write_text("sample,condition\n")

        # Create second config without metadata
        (tmp_path / "config_rna.yaml").write_text("test: config")

        config_files = find_assay_config_paths(tmp_path)
        metadata_files = find_metadata_paths(tmp_path)

        with pytest.raises(FileNotFoundError) as exc_info:
            validate_config_and_metadata(config_files, metadata_files)

        assert "Missing metadata file" in str(exc_info.value)
        assert "RNA" in str(exc_info.value)

    def test_skips_multiomics_config(self, tmp_path: Path):
        """Test that config_multiomics.yaml is skipped and doesn't require metadata."""
        # Create regular assay configs
        (tmp_path / "config_atac.yaml").write_text("test: config")
        (tmp_path / "metadata_atac.csv").write_text("sample,condition\n")

        # Create multiomics config WITHOUT metadata (should not raise error)
        (tmp_path / "config_multiomics.yaml").write_text("test: multiomics config")

        config_files = find_assay_config_paths(tmp_path)
        metadata_files = find_metadata_paths(tmp_path)

        # Should only find atac, not multiomics
        assert len(config_files) == 1
        assert Assay.ATAC in config_files
        assert Assay.MULTIOMICS not in config_files
        assert Assay.MULTIOMICS not in metadata_files


class TestMultiomicsConfig:
    """Tests for the MultiomicsConfig model."""

    def test_default_output_dir(self):
        """Test MultiomicsConfig with default output directory."""
        config = MultiomicsConfig()

        assert config.output_dir == "seqnado_output/"

    def test_custom_output_dir(self):
        """Test MultiomicsConfig with custom output directory."""
        config = MultiomicsConfig(output_dir="custom_output/")

        assert config.output_dir == "custom_output/"

    def test_default_binsize(self):
        """Test MultiomicsConfig with default binsize."""
        config = MultiomicsConfig()
        assert config.binsize == 1000

    def test_custom_binsize(self):
        """Test MultiomicsConfig with custom binsize."""
        config = MultiomicsConfig(binsize=500)
        assert config.binsize == 500

    def test_model_serialization(self):
        """Test that model can be serialized to dict."""
        config = MultiomicsConfig(output_dir="test_output/")

        data = config.model_dump()

        assert "output_dir" in data
        assert data["output_dir"] == "test_output/"
        assert data["binsize"] == 1000

    def test_model_from_dict(self):
        """Test that model can be created from dict."""
        data = {"output_dir": "custom/", "binsize": 200}
        config = MultiomicsConfig(**data)

        assert config.output_dir == "custom/"
        assert config.binsize == 200
