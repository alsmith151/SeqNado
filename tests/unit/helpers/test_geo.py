"""Tests for seqnado.workflow.helpers.geo module."""

from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from seqnado.workflow.helpers.geo import get_files_for_symlink, get_symlinked_files


class TestGetFilesForSymlink:
    """Tests for get_files_for_symlink function."""

    @patch("seqnado.outputs.files.GEOFiles")
    def test_returns_fastq_and_processed_files(self, mock_geo_files_class):
        """Test that function returns both fastq and processed files."""
        # Setup mocks
        output = Mock()
        output.files = [
            Path("/output/aligned/sample1.bam"),
            Path("/output/peaks/sample1.bed"),
        ]
        output.assay = "ChIP-seq"
        output.design_dataframe = pd.DataFrame()
        output.sample_names = ["sample1"]
        output.config = Mock()

        # Mock GEOFiles instance
        mock_geo_instance = Mock()
        mock_geo_instance.raw_files = {
            "sample1": ["sample1_R1.fastq.gz", "sample1_R2.fastq.gz"]
        }
        mock_geo_instance.processed_data_files = pd.DataFrame({
            "path": ["/output/peaks/sample1.bed"]
        })
        mock_geo_files_class.return_value = mock_geo_instance

        output_dir = "/output"

        result = get_files_for_symlink(output, output_dir)

        # Check fastq files are included
        assert "/output/fastqs/sample1_R1.fastq.gz" in result
        assert "/output/fastqs/sample1_R2.fastq.gz" in result
        # Check processed file is included
        assert "/output/peaks/sample1.bed" in result

    @patch("seqnado.outputs.files.GEOFiles")
    def test_excludes_geo_submission_files(self, mock_geo_files_class):
        """Test that geo_submission files are excluded from source files."""
        output = Mock()
        output.files = [
            Path("/output/aligned/sample1.bam"),
            Path("/output/geo_submission/sample1_symlink.bam"),  # Should be excluded
        ]
        output.assay = "ChIP-seq"
        output.design_dataframe = pd.DataFrame()
        output.sample_names = ["sample1"]
        output.config = Mock()

        mock_geo_instance = Mock()
        mock_geo_instance.raw_files = {"sample1": ["sample1_R1.fastq.gz"]}
        mock_geo_instance.processed_data_files = pd.DataFrame({"path": []})
        mock_geo_files_class.return_value = mock_geo_instance

        output_dir = "/output"

        # Call the function
        get_files_for_symlink(output, output_dir)

        # Verify GEOFiles was called with filtered source_files
        call_args = mock_geo_files_class.call_args
        source_files = call_args.kwargs["processed_files"]

        assert "/output/aligned/sample1.bam" in source_files
        assert "/output/geo_submission/sample1_symlink.bam" not in source_files

    @patch("seqnado.outputs.files.GEOFiles")
    def test_handles_empty_processed_files(self, mock_geo_files_class):
        """Test handling when there are no processed data files."""
        output = Mock()
        output.files = []
        output.assay = "ChIP-seq"
        output.design_dataframe = pd.DataFrame()
        output.sample_names = ["sample1"]
        output.config = Mock()

        mock_geo_instance = Mock()
        mock_geo_instance.raw_files = {"sample1": ["sample1_R1.fastq.gz"]}
        mock_geo_instance.processed_data_files = pd.DataFrame()  # Empty dataframe
        mock_geo_files_class.return_value = mock_geo_instance

        output_dir = "/output"

        result = get_files_for_symlink(output, output_dir)

        # Should only contain fastq files
        assert len(result) == 1
        assert "/output/fastqs/sample1_R1.fastq.gz" in result

    @patch("seqnado.outputs.files.GEOFiles")
    def test_sorted_fastq_files(self, mock_geo_files_class):
        """Test that fastq files are sorted."""
        output = Mock()
        output.files = []
        output.assay = "RNA-seq"
        output.design_dataframe = pd.DataFrame()
        output.sample_names = ["sample1", "sample2"]
        output.config = Mock()

        mock_geo_instance = Mock()
        mock_geo_instance.raw_files = {
            "sample2": ["sample2_R1.fastq.gz"],
            "sample1": ["sample1_R1.fastq.gz"],
        }
        mock_geo_instance.processed_data_files = pd.DataFrame()
        mock_geo_files_class.return_value = mock_geo_instance

        output_dir = "/output"

        result = get_files_for_symlink(output, output_dir)

        # Files should be sorted
        assert result == sorted(result)


class TestGetSymlinkedFiles:
    """Tests for get_symlinked_files function."""

    @patch("seqnado.outputs.files.GEOFiles")
    def test_returns_symlinked_paths(self, mock_geo_files_class):
        """Test that function returns paths in geo_submission directory."""
        output = Mock()
        output.files = [Path("/output/aligned/sample1.bam")]
        output.assay = "ChIP-seq"
        output.design_dataframe = pd.DataFrame()
        output.sample_names = ["sample1"]
        output.config = Mock()

        mock_geo_instance = Mock()
        mock_geo_instance.raw_files = {"sample1": ["sample1_R1.fastq.gz"]}
        mock_geo_instance.processed_data_files = pd.DataFrame({
            "output_file_name": ["sample1.bed"]
        })
        mock_geo_files_class.return_value = mock_geo_instance

        output_dir = "/output"

        result = get_symlinked_files(output, output_dir)

        # Check symlinked paths
        assert "/output/geo_submission/sample1_R1.fastq.gz" in result
        assert "/output/geo_submission/sample1.bed" in result

    @patch("seqnado.outputs.files.GEOFiles")
    def test_handles_multiple_fastq_pairs(self, mock_geo_files_class):
        """Test handling of paired-end fastq files."""
        output = Mock()
        output.files = []
        output.assay = "RNA-seq"
        output.design_dataframe = pd.DataFrame()
        output.sample_names = ["sample1"]
        output.config = Mock()

        mock_geo_instance = Mock()
        mock_geo_instance.raw_files = {
            "sample1": ["sample1_R1.fastq.gz", "sample1_R2.fastq.gz"]
        }
        mock_geo_instance.processed_data_files = pd.DataFrame()
        mock_geo_files_class.return_value = mock_geo_instance

        output_dir = "/output"

        result = get_symlinked_files(output, output_dir)

        assert len(result) == 2
        assert "/output/geo_submission/sample1_R1.fastq.gz" in result
        assert "/output/geo_submission/sample1_R2.fastq.gz" in result

    @patch("seqnado.outputs.files.GEOFiles")
    def test_empty_processed_files_dataframe(self, mock_geo_files_class):
        """Test handling of empty processed files dataframe."""
        output = Mock()
        output.files = []
        output.assay = "ChIP-seq"
        output.design_dataframe = pd.DataFrame()
        output.sample_names = ["sample1"]
        output.config = Mock()

        mock_geo_instance = Mock()
        mock_geo_instance.raw_files = {"sample1": ["sample1.fastq.gz"]}
        mock_geo_instance.processed_data_files = pd.DataFrame()  # Empty
        mock_geo_files_class.return_value = mock_geo_instance

        output_dir = "/output"

        result = get_symlinked_files(output, output_dir)

        # Should only have fastq file
        assert len(result) == 1
        assert result[0].endswith("sample1.fastq.gz")
