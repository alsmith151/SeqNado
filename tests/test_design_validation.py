"""Tests for the design.validation module."""
import pytest
import pandas as pd
from seqnado.design.validation import DesignDataFrame


@pytest.mark.unit
class TestDesignDataFrame:
    """Test cases for the DesignDataFrame validation schema."""

    def test_valid_dataframe_validation(self):
        """Test validation of a valid design DataFrame."""
        valid_data = pd.DataFrame({
            "sample_name": ["sample1", "sample2"],
            "r1": ["/path/to/sample1_R1.fastq.gz", "/path/to/sample2_R1.fastq.gz"],
            "r2": ["/path/to/sample1_R2.fastq.gz", "/path/to/sample2_R2.fastq.gz"],
            "norm_group": ["control", "treatment"]
        })
        
        # This should not raise an exception
        validated_df = DesignDataFrame.validate(valid_data)
        assert isinstance(validated_df, pd.DataFrame)
        assert len(validated_df) == 2

    def test_missing_required_columns(self):
        """Test validation fails with missing required columns."""
        invalid_data = pd.DataFrame({
            "sample_name": ["sample1"],
            # Missing r1 column
            "r2": ["/path/to/sample1_R2.fastq.gz"]
        })
        
        with pytest.raises(Exception):  # Should raise validation error
            DesignDataFrame.validate(invalid_data)

    def test_single_end_data(self):
        """Test validation of single-end data (no R2)."""
        single_end_data = pd.DataFrame({
            "sample_name": ["sample1", "sample2"],
            "r1": ["/path/to/sample1.fastq.gz", "/path/to/sample2.fastq.gz"],
            "r2": [None, None]  # No R2 reads
        })
        
        # This should not raise an exception
        validated_df = DesignDataFrame.validate(single_end_data)
        assert isinstance(validated_df, pd.DataFrame)

    @pytest.mark.parametrize("invalid_sample_name", [
        "",  # Empty string
        None,  # None value
        123,  # Non-string
    ])
    def test_invalid_sample_names(self, invalid_sample_name):
        """Test validation fails with invalid sample names."""
        invalid_data = pd.DataFrame({
            "sample_name": [invalid_sample_name],
            "r1": ["/path/to/file.fastq.gz"],
        })
        
        with pytest.raises(Exception):
            DesignDataFrame.validate(invalid_data)

    def test_duplicate_sample_names(self):
        """Test validation fails with duplicate sample names."""
        duplicate_data = pd.DataFrame({
            "sample_name": ["sample1", "sample1"],  # Duplicate
            "r1": ["/path/to/file1.fastq.gz", "/path/to/file2.fastq.gz"],
        })
        
        with pytest.raises(Exception):
            DesignDataFrame.validate(duplicate_data)

    def test_optional_metadata_columns(self):
        """Test validation with optional metadata columns."""
        data_with_metadata = pd.DataFrame({
            "sample_name": ["sample1", "sample2"],
            "r1": ["/path/to/sample1.fastq.gz", "/path/to/sample2.fastq.gz"],
            "norm_group": ["control", "treatment"],
            "condition": ["baseline", "drug_treatment"],
            "spikein_ratio": [None, 0.25]
        })
        
        validated_df = DesignDataFrame.validate(data_with_metadata)
        assert "condition" in validated_df.columns
        assert "spikein_ratio" in validated_df.columns

    @pytest.mark.edge_case
    def test_empty_dataframe(self):
        """Test validation of empty DataFrame."""
        empty_df = pd.DataFrame()
        
        with pytest.raises(Exception):
            DesignDataFrame.validate(empty_df)

    @pytest.mark.edge_case
    def test_very_long_paths(self):
        """Test validation with very long file paths."""
        long_path = "/very/long/path/" + "directory/" * 50 + "file.fastq.gz"
        data = pd.DataFrame({
            "sample_name": ["sample1"],
            "r1": [long_path],
        })
        
        # Should handle long paths without issues
        validated_df = DesignDataFrame.validate(data)
        assert len(validated_df) == 1
