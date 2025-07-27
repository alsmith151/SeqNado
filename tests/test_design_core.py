"""Tests for the design.core module."""
import pytest
from pydantic import ValidationError
from seqnado.inputs.core import Assay, PileupMethod, ScaleMethod, Metadata, clean_sample_name, extract_read_number, is_control_sample


# Test fixtures
@pytest.fixture
def sample_metadata():
    """Provide sample metadata for testing."""
    return {
        "norm_group": "treatment",
        "consensus_group": "high_dose",
        "deseq2": "condition"
    }


@pytest.fixture
def default_metadata():
    """Provide default metadata instance."""
    return Metadata()


@pytest.fixture
def custom_metadata(sample_metadata):
    """Provide custom metadata instance."""
    return Metadata(**sample_metadata)


@pytest.mark.unit
class TestAssay:
    """Test cases for the Assay enum."""

    def test_all_assay_values(self):
        """Test that all expected assay types are available."""
        expected_assays = {
            "rna", "atac", "snp", "chip", "cat", "meth", "mcc", "crispr"
        }
        actual_assays = {assay.value for assay in Assay}
        assert actual_assays == expected_assays

    def test_ip_assays(self):
        """Test identification of IP-requiring assays."""
        ip_assays = Assay.ip_assays()
        assert Assay.CHIP in ip_assays
        assert Assay.CAT in ip_assays
        assert len(ip_assays) == 2
        
        # Non-IP assays should not be in IP list
        assert Assay.RNA not in ip_assays
        assert Assay.ATAC not in ip_assays

    def test_non_ip_assays(self):
        """Test identification of non-IP assays."""
        non_ip_assays = Assay.non_ip_assays()
        assert Assay.RNA in non_ip_assays
        assert Assay.ATAC in non_ip_assays
        assert Assay.CHIP not in non_ip_assays
        assert Assay.CAT not in non_ip_assays

    @pytest.mark.parametrize("assay,expected_value", [
        (Assay.RNA, "rna"),
        (Assay.CHIP, "chip"),
        (Assay.ATAC, "atac"),
        (Assay.CAT, "cat"),
        (Assay.METH, "meth"),
        (Assay.MCC, "mcc"),
        (Assay.SNP, "snp"),
        (Assay.CRISPR, "crispr"),
    ])
    def test_enum_string_values(self, assay, expected_value):
        """Test that enums have correct string values."""
        assert assay.value == expected_value

    def test_enum_string_representation(self):
        """Test enum string representation."""
        assert str(Assay.ATAC) == "Assay.ATAC"
        assert repr(Assay.RNA) == "<Assay.RNA: 'rna'>"

    @pytest.mark.parametrize("assay", [
        Assay.CHIP,
        Assay.CAT,
    ])
    def test_ip_assay_membership(self, assay):
        """Test that IP assays are correctly identified."""
        assert assay in Assay.ip_assays()

    @pytest.mark.parametrize("assay", [
        Assay.RNA,
        Assay.ATAC,
        Assay.SNP,
        Assay.METH,
        Assay.MCC,
        Assay.CRISPR,
    ])
    def test_non_ip_assay_membership(self, assay):
        """Test that non-IP assays are correctly identified."""
        assert assay in Assay.non_ip_assays()


@pytest.mark.unit
class TestPileupMethod:
    """Test cases for the PileupMethod enum."""

    def test_pileup_methods(self):
        """Test that all expected pileup methods are available."""
        expected_methods = {"deeptools", "homer", "bamnado"}
        actual_methods = {method.value for method in PileupMethod}
        assert actual_methods == expected_methods

    @pytest.mark.parametrize("method,expected_value", [
        (PileupMethod.DEEPTOOLS, "deeptools"),
        (PileupMethod.HOMER, "homer"),
        (PileupMethod.BAMNADO, "bamnado"),
    ])
    def test_method_values(self, method, expected_value):
        """Test that pileup methods have correct values."""
        assert method.value == expected_value


@pytest.mark.unit
class TestScaleMethod:
    """Test cases for the ScaleMethod enum."""

    def test_scale_methods(self):
        """Test that all expected scale methods are available."""
        expected_methods = {"unscaled", "csaw", "cpm", "rpkm", "spikein", "merged"}
        actual_methods = {method.value for method in ScaleMethod}
        assert actual_methods == expected_methods

    @pytest.mark.parametrize("method,expected_value", [
        (ScaleMethod.UNSCALED, "unscaled"),
        (ScaleMethod.CSAW, "csaw"),
        (ScaleMethod.CPM, "cpm"),
        (ScaleMethod.RPKM, "rpkm"),
        (ScaleMethod.SPIKEIN, "spikein"),
        (ScaleMethod.MERGED, "merged"),
    ])
    def test_method_values(self, method, expected_value):
        """Test that scale methods have correct values."""
        assert method.value == expected_value


@pytest.mark.unit
class TestMetadata:
    """Test cases for the Metadata class."""

    def test_metadata_creation_with_defaults(self, default_metadata):
        """Test metadata creation with default values."""
        assert default_metadata.norm_group == "all"
        assert default_metadata.consensus_group is None
        assert default_metadata.deseq2 is None
        assert default_metadata.assay is None

    def test_metadata_creation_with_values(self, custom_metadata):
        """Test metadata creation with custom values."""
        assert custom_metadata.norm_group == "treatment"
        assert custom_metadata.consensus_group == "high_dose"
        assert custom_metadata.deseq2 == "condition"

    def test_metadata_dict_conversion(self, custom_metadata):
        """Test metadata conversion to dictionary."""
        data_dict = custom_metadata.model_dump()
        
        expected_keys = {"norm_group", "consensus_group", "deseq2", "assay"}
        assert set(data_dict.keys()) == expected_keys
        assert data_dict["norm_group"] == "treatment"
        assert data_dict["consensus_group"] == "high_dose"
        assert data_dict["deseq2"] == "condition"

    def test_metadata_exclude_none(self, default_metadata):
        """Test metadata dict conversion excluding None values."""
        data_dict = default_metadata.model_dump(exclude_none=True)
        
        assert "norm_group" in data_dict
        assert data_dict["norm_group"] == "all"
        assert "consensus_group" not in data_dict
        assert "deseq2" not in data_dict
        assert "assay" not in data_dict

    @pytest.mark.parametrize("norm_group,consensus_group,deseq2", [
        ("treatment", "low_dose", "condition"),
        ("treatment", "high_dose", "time_point"),
    ])
    def test_metadata_parametrized_creation(self, norm_group, consensus_group, deseq2):
        """Test metadata creation with valid parameter combinations."""
        metadata = Metadata(
            norm_group=norm_group,
            consensus_group=consensus_group,
            deseq2=deseq2
        )
        assert metadata.norm_group == norm_group
        assert metadata.consensus_group == consensus_group
        assert metadata.deseq2 == deseq2

    def test_metadata_equality(self):
        """Test metadata equality comparison."""
        metadata1 = Metadata(norm_group="test", consensus_group="A")
        metadata2 = Metadata(norm_group="test", consensus_group="A")
        metadata3 = Metadata(norm_group="test", consensus_group="B")
        
        assert metadata1 == metadata2
        assert metadata1 != metadata3

    def test_metadata_copy(self, custom_metadata):
        """Test metadata copying."""
        copied = custom_metadata.model_copy()
        assert copied == custom_metadata
        assert copied is not custom_metadata

    def test_metadata_with_assay(self):
        """Test metadata creation with assay field."""
        from seqnado.inputs.core import Assay
        metadata = Metadata(assay=Assay.RNA, norm_group="test")
        assert metadata.assay == Assay.RNA
        assert metadata.norm_group == "test"


@pytest.mark.unit
class TestUtilityFunctions:
    """Test cases for utility functions in core module."""

    @pytest.mark.parametrize("input_name,expected", [
        ("sample_1_R1_001", "sample_1"),
        ("test-sample_R2", "test-sample"),
        ("sample.fastq", "sample.fastq"),  # No matching pattern
        ("sample_001.fq.gz", "sample_001.fq.gz"),  # No matching pattern
        ("complex_sample_name_R1_001.fastq.gz", "complex_sample_name_R1_001.fastq.gz"),  # No extension removal
        ("simple", "simple"),
        ("", ""),
        ("sample_S1_L001_R1_001", "sample"),  # Multiple pattern matches
    ])
    def test_clean_sample_name(self, input_name, expected):
        """Test sample name cleaning function."""
        assert clean_sample_name(input_name) == expected

    @pytest.mark.parametrize("filename,expected", [
        ("sample_R1.fastq", 1),
        ("sample_R2.fastq", 2),
        ("sample_1.fastq", 1),
        ("sample_2.fastq", 2),
        ("sample.fastq", None),  # No read number pattern
        ("sample_R1_001.fastq.gz", 1),
        ("sample_R2_001.fq.gz", 2),
        ("no_read_number.fq", None),  # No read number pattern
    ])
    def test_extract_read_number(self, filename, expected):
        """Test read number extraction from filenames."""
        assert extract_read_number(filename) == expected

    @pytest.mark.parametrize("sample_name,expected", [
        ("input_sample", True),
        ("control_sample", True),
        ("igg_sample", True),
        ("mock_sample", True),  # Added mock
        ("sample_Input", True),  # Case insensitive
        ("sample_CONTROL", True),  # Case insensitive
        ("Sample_IgG", True),  # Case insensitive
        ("treatment_sample", False),
        ("chip_sample", False),
        ("experiment_sample", False),
        ("", False),  # Empty string
    ])
    def test_is_control_sample_comprehensive(self, sample_name, expected):
        """Comprehensive test for control sample detection."""
        assert is_control_sample(sample_name) == expected

    def test_is_control_sample_basic(self):
        """Test basic control sample detection."""
        assert is_control_sample("input_sample")
        assert is_control_sample("control_sample")
        assert is_control_sample("igg_sample")
        assert not is_control_sample("chip_sample")
        assert not is_control_sample("treatment_sample")

    @pytest.mark.parametrize("sample_name,expected", [
        ("sample_input", True),
        ("sample_control", True),
        ("sample_igg", True),
        ("sample_Input", True),  # Case insensitive
        ("sample_CONTROL", True),  # Case insensitive
        ("treatment_sample", False),
        ("chip_sample", False),
    ])
    def test_is_control_sample_parametrized(self, sample_name, expected):
        """Parametrized test for control sample detection."""
        assert is_control_sample(sample_name) == expected

    # Edge case tests
    @pytest.mark.edge_case
    class TestEdgeCases:
        """Test edge cases for utility functions."""

        def test_clean_sample_name_edge_cases(self):
            """Test edge cases for sample name cleaning."""
            assert clean_sample_name("") == ""
            assert clean_sample_name("R1") == "R1"  # Just read indicator
            assert clean_sample_name("_R1_") == "_R1"  # Partial match

        def test_extract_read_number_edge_cases(self):
            """Test edge cases for read number extraction."""
            assert extract_read_number("") is None  # No match for empty
            assert extract_read_number("_R1") == 1
            assert extract_read_number("R2") == None
            assert extract_read_number("no_numbers_here") is None  # No match

        def test_is_control_sample_edge_cases(self):
            """Test edge cases for control sample detection."""
            assert not is_control_sample("")  # Empty string
            assert is_control_sample("INPUT")  # All caps, matches input
            assert is_control_sample("prefix_input_suffix")  # Input in middle
