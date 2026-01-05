"""Unit tests for seqnado/outputs/multiomics.py helper functions."""

from pathlib import Path
import pytest

from seqnado import Assay
from seqnado.outputs.multiomics import (
    none_str_to_none,
    find_assay_config_paths,
    find_metadata_paths,
    validate_config_and_metadata,
)


class TestNoneStrToNone:
    """Tests for none_str_to_none validator function."""

    def test_none_string_lowercase(self):
        """Test that 'none' string is converted to None."""
        assert none_str_to_none("none") is None

    def test_none_string_uppercase(self):
        """Test that 'NONE' string is converted to None."""
        assert none_str_to_none("NONE") is None

    def test_none_string_mixed_case(self):
        """Test that 'NoNe' string is converted to None."""
        assert none_str_to_none("NoNe") is None

    def test_none_string_with_whitespace(self):
        """Test that '  none  ' string is converted to None."""
        assert none_str_to_none("  none  ") is None

    def test_non_none_string_unchanged(self):
        """Test that regular strings are unchanged."""
        assert none_str_to_none("hello") == "hello"

    def test_empty_string_unchanged(self):
        """Test that empty string is unchanged."""
        assert none_str_to_none("") == ""

    def test_non_string_unchanged(self):
        """Test that non-strings are unchanged."""
        assert none_str_to_none(42) == 42
        assert none_str_to_none(None) is None
        assert none_str_to_none([1, 2, 3]) == [1, 2, 3]


class TestFindAssayConfigPaths:
    """Tests for find_assay_config_paths function."""

    def test_finds_single_assay_config(self, tmp_path):
        """Test finding a single assay config file."""
        config_file = tmp_path / "config_atac.yaml"
        config_file.touch()

        result = find_assay_config_paths(tmp_path)

        assert Assay.ATAC in result
        assert result[Assay.ATAC] == config_file

    def test_finds_multiple_assay_configs(self, tmp_path):
        """Test finding multiple assay config files."""
        atac_config = tmp_path / "config_atac.yaml"
        rna_config = tmp_path / "config_rna.yaml"
        chip_config = tmp_path / "config_chip.yaml"

        atac_config.touch()
        rna_config.touch()
        chip_config.touch()

        result = find_assay_config_paths(tmp_path)

        assert len(result) == 3
        assert Assay.ATAC in result
        assert Assay.RNA in result
        assert Assay.CHIP in result
        assert result[Assay.ATAC] == atac_config
        assert result[Assay.RNA] == rna_config
        assert result[Assay.CHIP] == chip_config

    def test_skips_multiomics_config(self, tmp_path):
        """Test that multiomics config is skipped."""
        atac_config = tmp_path / "config_atac.yaml"
        multiomics_config = tmp_path / "config_multiomics.yaml"

        atac_config.touch()
        multiomics_config.touch()

        result = find_assay_config_paths(tmp_path)

        assert len(result) == 1
        assert Assay.ATAC in result
        assert Assay.MULTIOMICS not in result

    def test_empty_directory(self, tmp_path):
        """Test that empty directory returns empty dict."""
        result = find_assay_config_paths(tmp_path)
        assert result == {}

    def test_ignores_non_yaml_files(self, tmp_path):
        """Test that non-yaml files are ignored."""
        (tmp_path / "config_atac.txt").touch()
        (tmp_path / "config_rna.json").touch()
        (tmp_path / "config_chip.yaml").touch()

        result = find_assay_config_paths(tmp_path)

        assert len(result) == 1
        assert Assay.CHIP in result

    def test_ignores_files_without_config_prefix(self, tmp_path):
        """Test that files without 'config_' prefix are ignored."""
        (tmp_path / "atac.yaml").touch()
        (tmp_path / "metadata_atac.yaml").touch()
        (tmp_path / "config_atac.yaml").touch()

        result = find_assay_config_paths(tmp_path)

        assert len(result) == 1
        assert Assay.ATAC in result

    def test_all_assay_types(self, tmp_path):
        """Test finding configs for all non-multiomics assay types."""
        assays_to_test = [
            Assay.ATAC, Assay.RNA, Assay.CHIP, Assay.CAT,
            Assay.SNP, Assay.METH, Assay.MCC, Assay.CRISPR
        ]

        for assay in assays_to_test:
            config_file = tmp_path / f"config_{assay.clean_name}.yaml"
            config_file.touch()

        result = find_assay_config_paths(tmp_path)

        assert len(result) == len(assays_to_test)
        for assay in assays_to_test:
            assert assay in result


class TestFindMetadataPaths:
    """Tests for find_metadata_paths function."""

    def test_finds_single_metadata_file(self, tmp_path):
        """Test finding a single metadata file."""
        metadata_file = tmp_path / "metadata_atac.csv"
        metadata_file.touch()

        result = find_metadata_paths(tmp_path)

        assert Assay.ATAC in result
        assert result[Assay.ATAC] == metadata_file

    def test_finds_multiple_metadata_files(self, tmp_path):
        """Test finding multiple metadata files."""
        atac_metadata = tmp_path / "metadata_atac.csv"
        rna_metadata = tmp_path / "metadata_rna.csv"
        chip_metadata = tmp_path / "metadata_chip.csv"

        atac_metadata.touch()
        rna_metadata.touch()
        chip_metadata.touch()

        result = find_metadata_paths(tmp_path)

        assert len(result) == 3
        assert Assay.ATAC in result
        assert Assay.RNA in result
        assert Assay.CHIP in result
        assert result[Assay.ATAC] == atac_metadata
        assert result[Assay.RNA] == rna_metadata
        assert result[Assay.CHIP] == chip_metadata

    def test_skips_multiomics_metadata(self, tmp_path):
        """Test that multiomics metadata is skipped."""
        atac_metadata = tmp_path / "metadata_atac.csv"
        multiomics_metadata = tmp_path / "metadata_multiomics.csv"

        atac_metadata.touch()
        multiomics_metadata.touch()

        result = find_metadata_paths(tmp_path)

        assert len(result) == 1
        assert Assay.ATAC in result
        assert Assay.MULTIOMICS not in result

    def test_empty_directory(self, tmp_path):
        """Test that empty directory returns empty dict."""
        result = find_metadata_paths(tmp_path)
        assert result == {}

    def test_ignores_non_csv_files(self, tmp_path):
        """Test that non-csv files are ignored."""
        (tmp_path / "metadata_atac.txt").touch()
        (tmp_path / "metadata_rna.xlsx").touch()
        (tmp_path / "metadata_chip.csv").touch()

        result = find_metadata_paths(tmp_path)

        assert len(result) == 1
        assert Assay.CHIP in result

    def test_ignores_files_without_metadata_prefix(self, tmp_path):
        """Test that files without 'metadata_' prefix are ignored."""
        (tmp_path / "atac.csv").touch()
        (tmp_path / "config_atac.csv").touch()
        (tmp_path / "metadata_atac.csv").touch()

        result = find_metadata_paths(tmp_path)

        assert len(result) == 1
        assert Assay.ATAC in result


class TestValidateConfigAndMetadata:
    """Tests for validate_config_and_metadata function."""

    def test_valid_matching_configs_and_metadata(self):
        """Test validation passes when configs and metadata match."""
        config_paths = {
            Assay.ATAC: Path("config_atac.yaml"),
            Assay.RNA: Path("config_rna.yaml"),
        }
        metadata_paths = {
            Assay.ATAC: Path("metadata_atac.csv"),
            Assay.RNA: Path("metadata_rna.csv"),
        }

        # Should not raise any exception
        validate_config_and_metadata(config_paths, metadata_paths)

    def test_missing_config_raises_error(self):
        """Test that missing config file raises FileNotFoundError."""
        config_paths = {
            Assay.ATAC: Path("config_atac.yaml"),
        }
        metadata_paths = {
            Assay.ATAC: Path("metadata_atac.csv"),
            Assay.RNA: Path("metadata_rna.csv"),  # RNA has metadata but no config
        }

        with pytest.raises(FileNotFoundError, match="Missing config file for assay RNA"):
            validate_config_and_metadata(config_paths, metadata_paths)

    def test_missing_metadata_raises_error(self):
        """Test that missing metadata file raises FileNotFoundError."""
        config_paths = {
            Assay.ATAC: Path("config_atac.yaml"),
            Assay.RNA: Path("config_rna.yaml"),  # RNA has config but no metadata
        }
        metadata_paths = {
            Assay.ATAC: Path("metadata_atac.csv"),
        }

        with pytest.raises(FileNotFoundError, match="Missing metadata file for assay RNA"):
            validate_config_and_metadata(config_paths, metadata_paths)

    def test_empty_configs_and_metadata(self):
        """Test validation passes with empty dictionaries."""
        config_paths = {}
        metadata_paths = {}

        # Should not raise any exception
        validate_config_and_metadata(config_paths, metadata_paths)

    def test_multiple_missing_files(self):
        """Test error message when first mismatch is found."""
        config_paths = {
            Assay.RNA: Path("config_rna.yaml"),
        }
        metadata_paths = {
            Assay.ATAC: Path("metadata_atac.csv"),
        }

        # Will raise for one of the assays with missing file
        # (The specific assay depends on set iteration order, so just check FileNotFoundError is raised)
        with pytest.raises(FileNotFoundError):
            validate_config_and_metadata(config_paths, metadata_paths)

    def test_error_message_includes_assay_name(self):
        """Test that error message includes the specific assay name."""
        config_paths = {}
        metadata_paths = {Assay.ATAC: Path("metadata_atac.csv")}

        with pytest.raises(FileNotFoundError) as exc_info:
            validate_config_and_metadata(config_paths, metadata_paths)

        assert "ATAC" in str(exc_info.value)

    def test_error_message_includes_helpful_command(self):
        """Test that error message includes helpful seqnado command."""
        config_paths = {}
        metadata_paths = {Assay.ATAC: Path("metadata_atac.csv")}

        with pytest.raises(FileNotFoundError) as exc_info:
            validate_config_and_metadata(config_paths, metadata_paths)

        assert "seqnado design" in str(exc_info.value)

    def test_validation_with_single_assay(self):
        """Test validation with a single assay."""
        config_paths = {Assay.ATAC: Path("config_atac.yaml")}
        metadata_paths = {Assay.ATAC: Path("metadata_atac.csv")}

        # Should not raise
        validate_config_and_metadata(config_paths, metadata_paths)
