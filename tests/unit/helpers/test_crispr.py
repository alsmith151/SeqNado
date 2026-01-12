"""Tests for seqnado.workflow.helpers.crispr module."""

import json
import tempfile
from pathlib import Path
from unittest.mock import Mock

import pytest

from seqnado.workflow.helpers.crispr import get_cutadapt_adapter_args


class TestGetCutadaptAdapterArgs:
    """Tests for get_cutadapt_adapter_args function."""

    def test_no_adapter_file_returns_base_options(self):
        """Test that missing adapter file returns base config options."""
        wildcards = Mock(sample="sample1")
        config = Mock()
        config.third_party_tools.cutadapt.command_line_arguments = "--quality-cutoff 20 --minimum-length 30"
        output_dir = "/nonexistent/path"

        result = get_cutadapt_adapter_args(wildcards, config, output_dir)
        assert result == "--quality-cutoff 20 --minimum-length 30"

    def test_adapter_file_with_both_adapters(self):
        """Test adapter file with both R1 and R2 adapters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wildcards = Mock(sample="sample1")
            config = Mock()
            config.third_party_tools.cutadapt.command_line_arguments = "--quality-cutoff 20"
            output_dir = tmpdir

            # Create resources directory and adapter file
            resources_dir = Path(tmpdir) / "resources"
            resources_dir.mkdir()
            adapter_data = {
                "adapter_r1": "AGATCGGAAGAGC",
                "adapter_r2": "AGATCGGAAGAGC"
            }
            adapter_file = resources_dir / "sample1_adapters.json"
            adapter_file.write_text(json.dumps(adapter_data))

            result = get_cutadapt_adapter_args(wildcards, config, output_dir)

            assert "--quality-cutoff 20" in result
            assert "-g 'AGATCGGAAGAGC'" in result
            assert "-G 'AGATCGGAAGAGC'" in result

    def test_adapter_file_with_only_r1(self):
        """Test adapter file with only R1 adapter."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wildcards = Mock(sample="sample1")
            config = Mock()
            config.third_party_tools.cutadapt.command_line_arguments = "--quality-cutoff 20"
            output_dir = tmpdir

            # Create resources directory and adapter file with only R1
            resources_dir = Path(tmpdir) / "resources"
            resources_dir.mkdir()
            adapter_data = {"adapter_r1": "AGATCGGAAGAGC"}
            adapter_file = resources_dir / "sample1_adapters.json"
            adapter_file.write_text(json.dumps(adapter_data))

            result = get_cutadapt_adapter_args(wildcards, config, output_dir)

            assert "-g 'AGATCGGAAGAGC'" in result
            assert "-G" not in result

    def test_removes_existing_adapter_flags(self):
        """Test that existing adapter flags are removed from base options."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wildcards = Mock(sample="sample1")
            config = Mock()
            # Include various adapter flags that should be removed
            config.third_party_tools.cutadapt.command_line_arguments = (
                "-g OLDADAPTER -a OLDADAPTER -G OLDADAPTER -A OLDADAPTER --quality-cutoff 20"
            )
            output_dir = tmpdir

            # Create resources directory and adapter file
            resources_dir = Path(tmpdir) / "resources"
            resources_dir.mkdir()
            adapter_data = {
                "adapter_r1": "NEWR1",
                "adapter_r2": "NEWR2"
            }
            adapter_file = resources_dir / "sample1_adapters.json"
            adapter_file.write_text(json.dumps(adapter_data))

            result = get_cutadapt_adapter_args(wildcards, config, output_dir)

            # Old adapters should be gone
            assert "OLDADAPTER" not in result
            # New adapters should be present
            assert "-g 'NEWR1'" in result
            assert "-G 'NEWR2'" in result
            # Other options preserved
            assert "--quality-cutoff 20" in result

    def test_removes_adapter_flags_without_space(self):
        """Test removal of adapter flags attached to sequences without space."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wildcards = Mock(sample="sample1")
            config = Mock()
            config.third_party_tools.cutadapt.command_line_arguments = (
                "-gOLDADAPTER --quality-cutoff 20"
            )
            output_dir = tmpdir

            resources_dir = Path(tmpdir) / "resources"
            resources_dir.mkdir()
            adapter_data = {"adapter_r1": "NEWR1"}
            adapter_file = resources_dir / "sample1_adapters.json"
            adapter_file.write_text(json.dumps(adapter_data))

            result = get_cutadapt_adapter_args(wildcards, config, output_dir)

            assert "OLDADAPTER" not in result
            assert "-g 'NEWR1'" in result

    def test_invalid_json_returns_base_options(self):
        """Test that invalid JSON returns base config options."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wildcards = Mock(sample="sample1")
            config = Mock()
            config.third_party_tools.cutadapt.command_line_arguments = "--quality-cutoff 20"
            output_dir = tmpdir

            # Create resources directory and invalid JSON file
            resources_dir = Path(tmpdir) / "resources"
            resources_dir.mkdir()
            adapter_file = resources_dir / "sample1_adapters.json"
            adapter_file.write_text("{ invalid json }")

            result = get_cutadapt_adapter_args(wildcards, config, output_dir)
            assert result == "--quality-cutoff 20"

    def test_empty_adapter_values(self):
        """Test adapter file with null/empty adapter values."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wildcards = Mock(sample="sample1")
            config = Mock()
            config.third_party_tools.cutadapt.command_line_arguments = "--quality-cutoff 20"
            output_dir = tmpdir

            # Create resources directory and adapter file with null values
            resources_dir = Path(tmpdir) / "resources"
            resources_dir.mkdir()
            adapter_data = {"adapter_r1": None, "adapter_r2": None}
            adapter_file = resources_dir / "sample1_adapters.json"
            adapter_file.write_text(json.dumps(adapter_data))

            result = get_cutadapt_adapter_args(wildcards, config, output_dir)

            # Should have base options but no adapter flags
            assert "--quality-cutoff 20" in result
            assert "-g" not in result
            assert "-G" not in result
