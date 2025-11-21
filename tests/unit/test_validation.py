"""Tests for seqnado.inputs.validation module."""

import pandas as pd
import pytest
from pandera.errors import SchemaError

from seqnado.inputs.validation import DesignDataFrame, ViewpointsFile


class TestViewpointsFile:
    """Tests for ViewpointsFile validation schema."""

    def test_valid_viewpoints_file(self):
        """Test that a valid viewpoints file passes validation."""
        data = {
            "Chromosome": ["chr1", "chr2", "chr3"],
            "Start": [100, 200, 300],
            "End": [200, 300, 400],
            "Name": ["viewpoint1", "viewpoint2", "viewpoint3"],
            "Strand": ["+", "-", "+"],
            "Score": [10.0, 20.0, 30.0],
        }
        df = pd.DataFrame(data)

        validated = ViewpointsFile.validate(df)
        assert validated is not None
        assert len(validated) == 3

    def test_viewpoints_file_without_optional_columns(self):
        """Test that viewpoints file works without optional Strand and Score."""
        data = {
            "Chromosome": ["chr1", "chr2"],
            "Start": [100, 200],
            "End": [200, 300],
            "Name": ["viewpoint1", "viewpoint2"],
        }
        df = pd.DataFrame(data)

        validated = ViewpointsFile.validate(df)
        assert validated is not None
        assert len(validated) == 2

    def test_viewpoint_names_with_underscores_valid(self):
        """Test that viewpoint names with underscores are valid."""
        data = {
            "Chromosome": ["chr1"],
            "Start": [100],
            "End": [200],
            "Name": ["my_viewpoint_1"],
        }
        df = pd.DataFrame(data)

        validated = ViewpointsFile.validate(df)
        assert validated is not None

    def test_viewpoint_names_with_spaces_invalid(self):
        """Test that viewpoint names with spaces are invalid."""
        data = {
            "Chromosome": ["chr1"],
            "Start": [100],
            "End": [200],
            "Name": ["my viewpoint"],
        }
        df = pd.DataFrame(data)

        with pytest.raises(SchemaError, match="check_viewpoint_names"):
            ViewpointsFile.validate(df)

    def test_viewpoint_names_with_special_chars_invalid(self):
        """Test that viewpoint names with special characters are invalid."""
        data = {
            "Chromosome": ["chr1"],
            "Start": [100],
            "End": [200],
            "Name": ["viewpoint@1"],
        }
        df = pd.DataFrame(data)

        with pytest.raises(SchemaError, match="check_viewpoint_names"):
            ViewpointsFile.validate(df)

    def test_coercion_of_types(self):
        """Test that types are properly coerced."""
        data = {
            "Chromosome": [1, 2],  # Will be coerced to string
            "Start": ["100", "200"],  # Will be coerced to int
            "End": ["200", "300"],
            "Name": ["vp1", "vp2"],
        }
        df = pd.DataFrame(data)

        validated = ViewpointsFile.validate(df)
        assert validated["Chromosome"].dtype == object
        assert validated["Start"].dtype == "int64"


class TestDesignDataFrame:
    """Tests for DesignDataFrame validation schema."""

    def test_valid_design_dataframe(self):
        """Test that a valid design dataframe passes validation."""
        data = {
            "sample_id": ["sample1", "sample2", "sample3"],
            "r1": ["read1.fq", "read2.fq", "read3.fq"],
        }
        df = pd.DataFrame(data)

        validated = DesignDataFrame.validate(df)
        assert validated is not None
        assert len(validated) == 3

    def test_design_with_optional_columns(self):
        """Test design dataframe with all optional columns."""
        data = {
            "sample_id": ["sample1", "sample2"],
            "scaling_group": ["batch1", "batch1"],
            "consensus_group": ["group1", "group2"],
            "deseq2": ["condition1", "condition2"],
            "r1": ["read1.fq", "read2.fq"],
            "r2": ["read1_r2.fq", "read2_r2.fq"],
            "r1_control": ["ctrl1.fq", "ctrl2.fq"],
            "r2_control": ["ctrl1_r2.fq", "ctrl2_r2.fq"],
            "ip": ["ip1", "ip2"],
            "control": ["ctrl1", "ctrl2"],
            "assay": ["ChIP", "ATAC"],
        }
        df = pd.DataFrame(data)

        validated = DesignDataFrame.validate(df)
        assert validated is not None
        assert len(validated) == 2

    def test_unique_sample_ids(self):
        """Test that sample IDs must be unique."""
        data = {
            "sample_id": ["sample1", "sample1"],  # Duplicate
            "r1": ["read1.fq", "read2.fq"],
        }
        df = pd.DataFrame(data)

        with pytest.raises(SchemaError, match="check_sample_name_is_unique"):
            DesignDataFrame.validate(df)

    def test_unique_sample_id_plus_ip_combination(self):
        """Test that sample_id + ip combination must be unique."""
        data = {
            "sample_id": ["sample1", "sample1"],  # Same sample_id
            "ip": ["ip1", "ip2"],  # Different IP makes them unique
            "r1": ["read1.fq", "read2.fq"],
        }
        df = pd.DataFrame(data)

        validated = DesignDataFrame.validate(df)
        assert validated is not None

    def test_duplicate_sample_id_plus_ip_fails(self):
        """Test that duplicate sample_id + ip combination fails."""
        data = {
            "sample_id": ["sample1", "sample1"],
            "ip": ["ip1", "ip1"],  # Same combination
            "r1": ["read1.fq", "read2.fq"],
        }
        df = pd.DataFrame(data)

        with pytest.raises(SchemaError, match="check_sample_name_is_unique"):
            DesignDataFrame.validate(df)

    def test_sample_id_with_underscores_valid(self):
        """Test that sample IDs with underscores and hyphens are valid."""
        data = {
            "sample_id": ["sample_1", "sample-2", "sample_3-test"],
            "r1": ["read1.fq", "read2.fq", "read3.fq"],
        }
        df = pd.DataFrame(data)

        validated = DesignDataFrame.validate(df)
        assert validated is not None

    def test_sample_id_with_spaces_invalid(self):
        """Test that sample IDs with spaces are invalid."""
        data = {
            "sample_id": ["sample 1"],
            "r1": ["read1.fq"],
        }
        df = pd.DataFrame(data)

        with pytest.raises(SchemaError, match="check_sample_id"):
            DesignDataFrame.validate(df)

    def test_sample_id_with_special_chars_invalid(self):
        """Test that sample IDs with special characters are invalid."""
        data = {
            "sample_id": ["sample@1"],
            "r1": ["read1.fq"],
        }
        df = pd.DataFrame(data)

        with pytest.raises(SchemaError, match="check_sample_id"):
            DesignDataFrame.validate(df)

    def test_scaling_group_optional(self):
        """Test that scaling_group column is optional."""
        data = {
            "sample_id": ["sample1"],
            "r1": ["read1.fq"],
        }
        df = pd.DataFrame(data)

        # Should validate successfully even without scaling_group
        validated = DesignDataFrame.validate(df)
        assert "sample_id" in validated.columns

    def test_consensus_group_optional(self):
        """Test that consensus_group column is optional."""
        data = {
            "sample_id": ["sample1"],
            "r1": ["read1.fq"],
        }
        df = pd.DataFrame(data)

        # Should validate successfully even without consensus_group
        validated = DesignDataFrame.validate(df)
        assert "sample_id" in validated.columns

    def test_nullable_deseq2_column(self):
        """Test that deseq2 column can be None."""
        data = {
            "sample_id": ["sample1", "sample2"],
            "r1": ["read1.fq", "read2.fq"],
            "deseq2": [None, "condition1"],
        }
        df = pd.DataFrame(data)

        validated = DesignDataFrame.validate(df)
        assert validated is not None
        assert pd.isna(validated.loc[0, "deseq2"])

    def test_assay_column_with_valid_assay(self):
        """Test that assay column accepts valid Assay enum values."""
        data = {
            "sample_id": ["sample1", "sample2"],
            "r1": ["read1.fq", "read2.fq"],
            "assay": ["ChIP", "ATAC"],
        }
        df = pd.DataFrame(data)

        validated = DesignDataFrame.validate(df)
        assert validated is not None
