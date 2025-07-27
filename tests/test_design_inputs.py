"""Tests for the design.inputs module."""
import pytest
import pathlib
from design.collections import SampleGroup, SampleGroups, SampleCollection, IPSampleCollection
from seqnado.inputs.core import Assay, Metadata
from seqnado.inputs.fastq import FastqFile, FastqSet, FastqSetIP, FastqFileIP
from seqnado.inputs.experiment import ExperimentIP


@pytest.fixture
def test_fastq_paths():
    """Provide test FASTQ file paths."""
    base_path = pathlib.Path(__file__).parent / "data" / "fastq"
    return {
        "rna_1": base_path / "rna_1.fastq.gz",
        "rna_2": base_path / "rna_2.fastq.gz",
        "atac_1": base_path / "atac_1.fastq.gz",
        "atac_2": base_path / "atac_2.fastq.gz",
        "chip_mll_1": base_path / "chip-rx_MLL_1.fastq.gz",
        "chip_mll_2": base_path / "chip-rx_MLL_2.fastq.gz",
        "chip_input_1": base_path / "chip-rx_input_1.fastq.gz",
        "chip_input_2": base_path / "chip-rx_input_2.fastq.gz",
    }


@pytest.fixture
def mock_design(test_fastq_paths):
    """Provide a real Design object for testing."""
    # Create FastqSet objects
    rna_fastq_set = FastqSet(
        sample_id="sample1",
        r1=FastqFile(path=test_fastq_paths["rna_1"]),
        r2=FastqFile(path=test_fastq_paths["rna_2"])
    )
    
    atac_fastq_set = FastqSet(
        sample_id="sample2", 
        r1=FastqFile(path=test_fastq_paths["atac_1"]),
        r2=FastqFile(path=test_fastq_paths["atac_2"])
    )
    
    # Single-end sample for variety
    single_fastq_set = FastqSet(
        sample_id="sample3",
        r1=FastqFile(path=test_fastq_paths["rna_1"]),
        r2=None
    )
    
    # Create metadata
    metadata = [
        Metadata(norm_group="all"),
        Metadata(norm_group="all"),
        Metadata(norm_group="all")
    ]
    
    return SampleCollection(
        assay=Assay.RNA,
        fastq_sets=[rna_fastq_set, atac_fastq_set, single_fastq_set],
        metadata=metadata
    )


@pytest.fixture
def mock_design_with_groups(test_fastq_paths):
    """Provide a real Design object with norm groups."""
    # Create FastqSet objects
    control_fastq_set = FastqSet(
        sample_id="sample1",
        r1=FastqFile(path=test_fastq_paths["rna_1"]),
        r2=FastqFile(path=test_fastq_paths["rna_2"])
    )
    
    treatment1_fastq_set = FastqSet(
        sample_id="sample2", 
        r1=FastqFile(path=test_fastq_paths["atac_1"]),
        r2=FastqFile(path=test_fastq_paths["atac_2"])
    )
    
    treatment2_fastq_set = FastqSet(
        sample_id="sample3",
        r1=FastqFile(path=test_fastq_paths["rna_1"]),
        r2=None
    )
    
    # Create metadata with different norm groups
    metadata = [
        Metadata(norm_group="control"),
        Metadata(norm_group="treatment"),
        Metadata(norm_group="treatment")
    ]
    
    return SampleCollection(
        assay=Assay.RNA,
        fastq_sets=[control_fastq_set, treatment1_fastq_set, treatment2_fastq_set],
        metadata=metadata
    )


@pytest.fixture
def mock_design_ip(test_fastq_paths):
    """Provide a real DesignIP object for testing."""
    # Create IP FastqSetIP objects
    ip_experiment1 = ExperimentIP(
        ip=FastqSetIP(
            sample_id="exp1",
            r1=FastqFileIP(path=test_fastq_paths["chip_mll_1"], ip="chip1", is_control=False),
            r2=FastqFileIP(path=test_fastq_paths["chip_mll_2"], ip="chip1", is_control=False),
            antibody="chip1"
        ),
        control=FastqSetIP(
            sample_id="exp1",
            r1=FastqFileIP(path=test_fastq_paths["chip_input_1"], ip="input", is_control=True),
            r2=FastqFileIP(path=test_fastq_paths["chip_input_2"], ip="input", is_control=True),
            antibody="input"
        )
    )
    
    ip_experiment2 = ExperimentIP(
        ip=FastqSetIP(
            sample_id="exp2",
            r1=FastqFileIP(path=test_fastq_paths["chip_mll_1"], ip="chip2", is_control=False),
            r2=FastqFileIP(path=test_fastq_paths["chip_mll_2"], ip="chip2", is_control=False),
            antibody="chip2"
        ),
        control=FastqSetIP(
            sample_id="exp2",
            r1=FastqFileIP(path=test_fastq_paths["chip_input_1"], ip="input", is_control=True),
            r2=FastqFileIP(path=test_fastq_paths["chip_input_2"], ip="input", is_control=True),
            antibody="input"
        )
    )
    
    # Create metadata
    metadata = [
        Metadata(norm_group="all"),
        Metadata(norm_group="all")
    ]
    
    return IPSampleCollection(
        assay=Assay.CHIP,
        experiments=[ip_experiment1, ip_experiment2],
        metadata=metadata
    )


@pytest.mark.unit
class TestSampleGroup:
    """Test cases for the SampleGroup class."""

    def test_sample_group_creation(self):
        """Test basic SampleGroup creation."""
        group = SampleGroup(
            group="treatment",
            samples=["sample1", "sample2"],
            reference_sample="sample1"
        )
        assert group.group == "treatment"
        assert group.samples == ["sample1", "sample2"]
        assert group.reference_sample == "sample1"

    def test_norm_group_defaults(self):
        """Test SampleGroup creation with defaults."""
        group = SampleGroup(samples=["sample1"])
        assert group.group == "all"
        assert group.reference_sample is None

    def test_norm_group_len(self):
        """Test SampleGroup length."""
        group = SampleGroup(samples=["sample1", "sample2", "sample3"])
        assert len(group) == 3

    def test_norm_group_contains(self):
        """Test SampleGroup membership testing."""
        group = SampleGroup(samples=["sample1", "sample2"])
        assert "sample1" in group
        assert "sample2" in group
        assert "sample3" not in group

    def test_norm_group_str(self):
        """Test SampleGroup string representation."""
        group = SampleGroup(
            group="test", 
            samples=["s1", "s2"], 
            reference_sample="s1"
        )
        expected = "SampleGroup(group='test', samples=2, reference='s1')"
        assert str(group) == expected

    def test_from_design_basic(self, mock_design):
        """Test creating SampleGroup from basic design."""
        group = SampleGroup.from_design(mock_design)
        assert group.group == "all"
        assert set(group.samples) == {"sample1", "sample2", "sample3"}
        assert group.reference_sample == "sample1"  # First sample

    def test_from_design_with_subset(self, mock_design_with_groups):
        """Test creating SampleGroup with subset values."""
        group = SampleGroup.from_design(
            mock_design_with_groups,
            subset_column="norm_group",
            subset_value=["treatment"]
        )
        assert group.group == "treatment"
        assert set(group.samples) == {"sample2", "sample3"}

    def test_from_design_ip_no_controls(self, mock_design_ip):
        """Test creating SampleGroup from DesignIP without controls."""
        group = SampleGroup.from_design(
            mock_design_ip,
            include_controls=False
        )
        # Should include IP samples with naming convention
        assert len(group.samples) > 0

    def test_from_design_ip_with_controls(self, mock_design_ip):
        """Test creating SampleGroup from DesignIP with controls."""
        group = SampleGroup.from_design(
            mock_design_ip,
            include_controls=True
        )
        # Should include both IP and control samples
        assert len(group.samples) > 0

    def test_from_design_empty_result_raises_error(self, mock_design):
        """Test that empty results raise ValueError."""
        with pytest.raises(ValueError, match="No samples found"):
            SampleGroup.from_design(
                mock_design,
                subset_column="nonexistent",
                subset_value=["nonexistent"]
            )

    def test_prepare_dataframe_design(self, mock_design):
        """Test _prepare_dataframe with Design object."""
        df = SampleGroup._prepare_dataframe(mock_design, include_controls=False)
        assert "sample_fullname" in df.index.names or df.index.name == "sample_fullname"
        assert len(df) == 3

    def test_prepare_dataframe_design_ip(self, mock_design_ip):
        """Test _prepare_dataframe with DesignIP object."""
        df = SampleGroup._prepare_dataframe(mock_design_ip, include_controls=False)
        assert len(df) >= 2  # Should have IP samples


@pytest.mark.unit
class TestSampleGroups:
    """Test cases for the SampleGroups class."""

    def test_norm_groups_creation(self):
        """Test basic SampleGroups creation."""
        group1 = SampleGroup(group="control", samples=["sample1"])
        group2 = SampleGroup(group="treatment", samples=["sample2", "sample3"])
        groups = SampleGroups(groups=[group1, group2])
        
        assert len(groups) == 2
        assert len(groups.groups) == 2

    def test_from_design_no_groups_column(self, mock_design):
        """Test creating SampleGroups when no grouping column exists."""
        groups = SampleGroups.from_sample_collection(mock_design)
        assert len(groups) == 1
        assert groups.groups[0].group == "all"
        assert set(groups.groups[0].samples) == {"sample1", "sample2", "sample3"}

    def test_from_design_with_groups_column(self, mock_design_with_groups):
        """Test creating SampleGroups with grouping column."""
        groups = SampleGroups.from_sample_collection(
            mock_design_with_groups,
            subset_column="norm_group"
        )
        assert len(groups) == 2  # control and treatment
        
        group_names = {group.group for group in groups.groups}
        assert group_names == {"control", "treatment"}

    def test_sample_to_group_mapping(self):
        """Test sample to group mapping property."""
        group1 = SampleGroup(group="A", samples=["sample1", "sample2"])
        group2 = SampleGroup(group="B", samples=["sample3"])
        groups = SampleGroups(groups=[group1, group2])
        
        mapping = groups.sample_to_group_mapping
        assert mapping["sample1"] == "A"
        assert mapping["sample2"] == "A"
        assert mapping["sample3"] == "B"

    def test_group_to_samples_mapping(self):
        """Test group to samples mapping property."""
        group1 = SampleGroup(group="A", samples=["sample1", "sample2"])
        group2 = SampleGroup(group="B", samples=["sample3"])
        groups = SampleGroups(groups=[group1, group2])
        
        mapping = groups.group_to_samples_mapping
        assert mapping["A"] == ["sample1", "sample2"]
        assert mapping["B"] == ["sample3"]

    def test_get_sample_group(self):
        """Test getting group for a sample."""
        group1 = SampleGroup(group="control", samples=["sample1"])
        group2 = SampleGroup(group="treatment", samples=["sample2"])
        groups = SampleGroups(groups=[group1, group2])
        
        assert groups.get_sample_group("sample1") == "control"
        assert groups.get_sample_group("sample2") == "treatment"

    def test_get_sample_group_not_found(self):
        """Test getting group for nonexistent sample raises KeyError."""
        group = SampleGroup(group="test", samples=["sample1"])
        groups = SampleGroups(groups=[group])
        
        with pytest.raises(KeyError, match="Sample 'nonexistent' not found"):
            groups.get_sample_group("nonexistent")

    def test_get_samples_in_group(self):
        """Test getting samples in a group."""
        group1 = SampleGroup(group="A", samples=["sample1", "sample2"])
        group2 = SampleGroup(group="B", samples=["sample3"])
        groups = SampleGroups(groups=[group1, group2])
        
        assert groups.get_samples_in_group("A") == ["sample1", "sample2"]
        assert groups.get_samples_in_group("B") == ["sample3"]

    def test_get_samples_in_group_not_found(self):
        """Test getting samples for nonexistent group raises KeyError."""
        group = SampleGroup(group="test", samples=["sample1"])
        groups = SampleGroups(groups=[group])
        
        with pytest.raises(KeyError, match="Group 'nonexistent' not found"):
            groups.get_samples_in_group("nonexistent")

    def test_get_group_by_name(self):
        """Test getting SampleGroup object by name."""
        group1 = SampleGroup(group="control", samples=["sample1"])
        group2 = SampleGroup(group="treatment", samples=["sample2"])
        groups = SampleGroups(groups=[group1, group2])
        
        retrieved = groups.get_group_by_name("control")
        assert retrieved == group1
        assert retrieved.samples == ["sample1"]

    def test_get_group_by_name_not_found(self):
        """Test getting nonexistent group by name raises KeyError."""
        group = SampleGroup(group="test", samples=["sample1"])
        groups = SampleGroups(groups=[group])
        
        with pytest.raises(KeyError, match="Group 'nonexistent' not found"):
            groups.get_group_by_name("nonexistent")

    def test_groups_iteration(self):
        """Test iterating over groups."""
        group1 = SampleGroup(group="A", samples=["sample1"])
        group2 = SampleGroup(group="B", samples=["sample2"])
        groups = SampleGroups(groups=[group1, group2])
        
        iterated = list(groups)
        assert len(iterated) == 2
        assert group1 in iterated
        assert group2 in iterated

    def test_groups_contains(self):
        """Test group membership testing."""
        group1 = SampleGroup(group="control", samples=["sample1"])
        group2 = SampleGroup(group="treatment", samples=["sample2"])
        groups = SampleGroups(groups=[group1, group2])
        
        assert "control" in groups
        assert "treatment" in groups
        assert "nonexistent" not in groups

    def test_groups_str(self):
        """Test SampleGroups string representation."""
        group1 = SampleGroup(group="A", samples=["s1", "s2"])
        group2 = SampleGroup(group="B", samples=["s3"])
        groups = SampleGroups(groups=[group1, group2])
        
        result = str(groups)
        assert "SampleGroups(2 groups:" in result
        assert "'A': 2 samples" in result
        assert "'B': 1 samples" in result

    # Test backward compatibility properties
    def test_backward_compatibility_sample_groups(self):
        """Test backward compatibility sample_groups property."""
        group1 = SampleGroup(group="A", samples=["sample1"])
        groups = SampleGroups(groups=[group1])
        
        # Should be alias for group_to_samples_mapping
        assert groups.sample_groups == groups.group_to_samples_mapping

    def test_backward_compatibility_group_samples(self):
        """Test backward compatibility group_samples property."""
        group1 = SampleGroup(group="A", samples=["sample1"])
        groups = SampleGroups(groups=[group1])
        
        # Should be alias for sample_to_group_mapping
        assert groups.group_samples == groups.sample_to_group_mapping

    def test_backward_compatibility_get_grouped_samples(self):
        """Test backward compatibility get_grouped_samples method."""
        group1 = SampleGroup(group="A", samples=["sample1", "sample2"])
        groups = SampleGroups(groups=[group1])
        
        # Should be alias for get_samples_in_group
        assert groups.get_grouped_samples("A") == groups.get_samples_in_group("A")


@pytest.mark.integration
class TestSampleGroupsIntegration:
    """Integration tests for SampleGroups with real design objects."""
    
    @pytest.fixture
    def sample_design_data(self, tmp_path):
        """Create sample design data for testing."""
        # This would require actual Design/DesignIP objects
        # For now, we'll skip this and focus on unit tests
        pytest.skip("Integration tests require full design objects")

    def test_integration_with_real_design(self, sample_design_data):
        """Test SampleGroups with real Design objects."""
        # This would test the full integration
        pass


@pytest.mark.parametrize("group_name,samples,reference", [
    ("control", ["c1", "c2"], "c1"),
    ("treatment", ["t1", "t2", "t3"], "t1"),
    (1, ["sample1"], "sample1"),
    ("mixed", ["a", "b", "c", "d"], "b"),
])
def test_norm_group_parametrized_creation(group_name, samples, reference):
    """Parametrized test for SampleGroup creation."""
    group = SampleGroup(
        group=group_name,
        samples=samples,
        reference_sample=reference
    )
    assert group.group == group_name
    assert group.samples == samples
    assert group.reference_sample == reference
    assert len(group) == len(samples)


@pytest.mark.edge_case
class TestSampleGroupsEdgeCases:
    """Test edge cases for SampleGroup and SampleGroups."""

    def test_empty_samples_list(self):
        """Test SampleGroup with empty samples list."""
        group = SampleGroup(samples=[])
        assert len(group) == 0
        assert "anything" not in group

    def test_single_sample_group(self):
        """Test SampleGroup with single sample."""
        group = SampleGroup(samples=["only_sample"])
        assert len(group) == 1
        assert "only_sample" in group

    def test_numeric_group_names(self):
        """Test SampleGroup with numeric group names."""
        group = SampleGroup(group=123, samples=["sample1"])
        assert group.group == 123
        assert isinstance(group.group, int)

    def test_empty_norm_groups(self):
        """Test SampleGroups with no groups."""
        groups = SampleGroups(groups=[])
        assert len(groups) == 0
        assert list(groups) == []
        assert groups.sample_to_group_mapping == {}
        assert groups.group_to_samples_mapping == {}

    def test_duplicate_samples_across_groups(self):
        """Test behavior with duplicate samples across groups."""
        group1 = SampleGroup(group="A", samples=["sample1", "sample2"])
        group2 = SampleGroup(group="B", samples=["sample1", "sample3"])  # sample1 duplicated
        groups = SampleGroups(groups=[group1, group2])
        
        # The mapping should contain the last occurrence
        mapping = groups.sample_to_group_mapping
        # This documents current behavior - last group wins
        assert mapping["sample1"] == "B"
        assert mapping["sample2"] == "A"
        assert mapping["sample3"] == "B"
