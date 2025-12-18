"""Unit tests for seqnado.config.multiomics module."""

from pathlib import Path

import pytest

from seqnado.outputs.multiomics import (
    find_assay_config_paths,
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

        config_files, metadata_files = find_assay_config_paths(tmp_path)

        assert len(config_files) == 1
        assert len(metadata_files) == 1
        assert "chip" in config_files
        assert "chip" in metadata_files
        assert config_files["chip"]["path"] == str(config_file)
        assert metadata_files["chip"] == str(metadata_file)

    def test_finds_multiple_assay_configs(self, tmp_path: Path):
        """Test finding multiple assay configs and metadata pairs."""
        assays = ["chip", "rna", "atac"]

        for assay in assays:
            config_file = tmp_path / f"config_{assay}.yaml"
            config_file.write_text(f"assay: {assay}")

            metadata_file = tmp_path / f"metadata_{assay}.csv"
            metadata_file.write_text("sample,condition\n")

        config_files, metadata_files = find_assay_config_paths(tmp_path)

        assert len(config_files) == 3
        assert len(metadata_files) == 3

        for assay in assays:
            assert assay in config_files
            assert assay in metadata_files
            assert config_files[assay]["path"] == str(tmp_path / f"config_{assay}.yaml")
            assert metadata_files[assay] == str(tmp_path / f"metadata_{assay}.csv")

    def test_raises_error_when_metadata_missing(self, tmp_path: Path):
        """Test that FileNotFoundError is raised when metadata file is missing."""
        # Create config file but not metadata
        config_file = tmp_path / "config_chip.yaml"
        config_file.write_text("test: config")

        with pytest.raises(FileNotFoundError) as exc_info:
            find_assay_config_paths(tmp_path)

        assert "Missing metadata file" in str(exc_info.value)
        assert "metadata_chip.csv" in str(exc_info.value)
        assert "seqnado design chip" in str(exc_info.value)

    def test_empty_directory_returns_empty_dicts(self, tmp_path: Path):
        """Test that empty directory returns empty dictionaries."""
        config_files, metadata_files = find_assay_config_paths(tmp_path)

        assert config_files == {}
        assert metadata_files == {}

    def test_ignores_non_config_yaml_files(self, tmp_path: Path):
        """Test that non-config YAML files are ignored."""
        # Create files that don't match the config_*.yaml pattern
        (tmp_path / "other.yaml").write_text("test: data")
        (tmp_path / "settings.yaml").write_text("test: data")
        (tmp_path / "config.yaml").write_text("test: data")  # Missing underscore

        config_files, metadata_files = find_assay_config_paths(tmp_path)

        assert config_files == {}
        assert metadata_files == {}

    def test_handles_assay_names_with_hyphens(self, tmp_path: Path):
        """Test handling assay names with hyphens."""
        config_file = tmp_path / "config_chip-rx.yaml"
        config_file.write_text("test: config")

        metadata_file = tmp_path / "metadata_chip-rx.csv"
        metadata_file.write_text("sample,condition\n")

        config_files, metadata_files = find_assay_config_paths(tmp_path)

        assert "chip-rx" in config_files
        assert "chip-rx" in metadata_files

    def test_partial_configs_raise_error(self, tmp_path: Path):
        """Test that having some configs with metadata and some without raises error."""
        # Create first pair correctly
        (tmp_path / "config_chip.yaml").write_text("test: config")
        (tmp_path / "metadata_chip.csv").write_text("sample,condition\n")

        # Create second config without metadata
        (tmp_path / "config_rna.yaml").write_text("test: config")

        with pytest.raises(FileNotFoundError) as exc_info:
            find_assay_config_paths(tmp_path)

        assert "metadata_rna.csv" in str(exc_info.value)

    def test_skips_multiomics_config(self, tmp_path: Path):
        """Test that config_multiomics.yaml is skipped and doesn't require metadata."""
        # Create regular assay configs
        (tmp_path / "config_atac.yaml").write_text("test: config")
        (tmp_path / "metadata_atac.csv").write_text("sample,condition\n")

        # Create multiomics config WITHOUT metadata (should not raise error)
        (tmp_path / "config_multiomics.yaml").write_text("test: multiomics config")

        config_files, metadata_files = find_assay_config_paths(tmp_path)

        # Should only find atac, not multiomics
        assert len(config_files) == 1
        assert "atac" in config_files
        assert "multiomics" not in config_files
        assert "multiomics" not in metadata_files


# class TestMultiomicsOutput:
#     """Tests for the MultiomicsOutput model."""

#     def test_default_output_dir(self):
#         """Test MultiomicsOutput with default output directory."""
#         output = MultiomicsConfig()

#         assert output.output_dir == "seqnado_output/"

#     def test_custom_output_dir(self):
#         """Test MultiomicsConfig with custom output directory."""
#         output = MultiomicsConfig(output_dir="custom_output/")

#         assert output.output_dir == "custom_output/"

#     def test_none_string_converts_to_none(self):
#         """Test that 'none' string is converted to None for output_dir."""
#         output = MultiomicsConfig(output_dir="none")

#         assert output.output_dir is None

#     def test_uppercase_none_string_converts_to_none(self):
#         """Test that 'NONE' string is converted to None for output_dir."""
#         output = MultiomicsConfig(output_dir="NONE")

#         assert output.output_dir is None

#     def test_summary_report_property(self):
#         """Test summary_report property returns correct path."""
#         output = MultiomicsConfig(output_dir="test_output/")

#         assert output.summary_report == str(
#             Path("test_output/") / "multiomics_summary.txt"
#         )

#     def test_heatmap_property(self):
#         """Test heatmap property returns correct path."""
#         output = MultiomicsConfig(output_dir="test_output/")

#         expected = str(Path("test_output/") / "multiomics" / "heatmap" / "heatmap.pdf")
#         assert output.heatmap == expected

#     def test_metaplot_property(self):
#         """Test metaplot property returns correct path."""
#         output = MultiomicsConfig(output_dir="test_output/")

#         expected = str(Path("test_output/") / "multiomics" / "heatmap" / "metaplot.pdf")
#         assert output.metaplot == expected

#     def test_all_outputs_property(self):
#         """Test all_outputs property returns list of all output files."""
#         output = MultiomicsConfig(output_dir="test_output/")

#         all_outputs = output.all_outputs

#         assert len(all_outputs) == 4
#         assert output.summary_report in all_outputs
#         assert output.heatmap in all_outputs
#         assert output.metaplot in all_outputs

#     def test_all_outputs_order(self):
#         """Test that all_outputs returns files in expected order."""
#         output = MultiomicsConfig(output_dir="test_output/")

#         all_outputs = output.all_outputs

#         assert all_outputs[0] == output.summary_report
#         assert all_outputs[1] == output.heatmap
#         assert all_outputs[2] == output.metaplot

#     def test_paths_with_trailing_slash(self):
#         """Test that paths work correctly with trailing slash."""
#         output = MultiomicsConfig(output_dir="test_output/")

#         assert "test_output/" in output.summary_report
#         assert "test_output/" in output.heatmap

#     def test_paths_without_trailing_slash(self):
#         """Test that paths work correctly without trailing slash."""
#         output = MultiomicsConfig(output_dir="test_output")

#         assert "test_output" in output.summary_report
#         assert "test_output" in output.heatmap

#     def test_relative_paths(self):
#         """Test that relative paths work correctly."""
#         output = MultiomicsConfig(output_dir="./results/multiomics/")

#         # Path normalizes ./results to results, so check for the path components
#         assert "results/multiomics" in output.summary_report
#         assert "multiomics_summary.txt" in output.summary_report

#     def test_absolute_paths(self):
#         """Test that absolute paths work correctly."""
#         output = MultiomicsConfig(output_dir="/absolute/path/output/")

#         assert output.summary_report.startswith("/absolute/path/output")
#         assert output.heatmap.startswith("/absolute/path/output")

#     def test_model_is_immutable_after_creation(self):
#         """Test that model fields can be accessed after creation."""
#         output = MultiomicsConfig(output_dir="test/")

#         # Properties should work multiple times
#         assert output.summary_report == output.summary_report
#         assert output.heatmap == output.heatmap
#         assert output.metaplot == output.metaplot
#         assert output.all_outputs == output.all_outputs

#     def test_model_serialization(self):
#         """Test that model can be serialized to dict."""
#         output = MultiomicsConfig(output_dir="test_output/")

#         data = output.model_dump()

#         assert "output_dir" in data
#         assert data["output_dir"] == "test_output/"

#     def test_model_from_dict(self):
#         """Test that model can be created from dict."""
#         data = {"output_dir": "custom/"}
#         output = MultiomicsConfig(**data)

#         assert output.output_dir == "custom/"

#     def test_none_output_dir_causes_error(self):
#         """Test that None output_dir causes TypeError when accessing path properties."""
#         output = MultiomicsConfig(output_dir=None)

#         # Path(None) raises TypeError, so accessing properties should fail
#         with pytest.raises(TypeError):
#             _ = output.summary_report

#         with pytest.raises(TypeError):
#             _ = output.heatmap

#         with pytest.raises(TypeError):
#             _ = output.metaplot
