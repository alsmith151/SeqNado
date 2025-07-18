"""Tests for the design.inputs module."""
import pytest
import pandas as pd
from seqnado.design.inputs import NormGroup, NormGroups


class MockDesign:
    """Mock Design class for testing."""
    
    def __init__(self, sample_names, metadata_dict=None):
        self.sample_names = sample_names
        self._metadata_dict = metadata_dict or {}
    
    def to_dataframe(self):
        """Return a mock DataFrame."""
        data = []
        for name in self.sample_names:
            row = {"sample_name": name}
            row.update(self._metadata_dict.get(name, {}))
            data.append(row)
        return pd.DataFrame(data)


class MockDesignIP:
    """Mock DesignIP class for testing."""
    
    def __init__(self, sample_names, ip_names=None, control_names=None):
        self.sample_names = sample_names
        self._ip_names = ip_names or []
        self._control_names = control_names or []
    
    def to_dataframe(self):
        """Return a mock DataFrame."""
        data = []
        for i, name in enumerate(self.sample_names):
            ip_name = self._ip_names[i] if i < len(self._ip_names) else "chip"
            control_name = self._control_names[i] if i < len(self._control_names) else "input"
            data.append({
                "sample_name": name,
                "ip": ip_name,
                "control": control_name
            })
        return pd.DataFrame(data)


@pytest.fixture
def mock_design():
    """Provide a mock Design object."""
    return MockDesign(["sample1", "sample2", "sample3"])


@pytest.fixture
def mock_design_with_groups():
    """Provide a mock Design object with norm groups."""
    metadata = {
        "sample1": {"norm_group": "control"},
        "sample2": {"norm_group": "treatment"},
        "sample3": {"norm_group": "treatment"}
    }
    design = MockDesign(["sample1", "sample2", "sample3"], metadata)
    
    # Add the norm_group column to the DataFrame
    def to_dataframe_with_groups():
        df = pd.DataFrame([
            {"sample_name": "sample1", "norm_group": "control"},
            {"sample_name": "sample2", "norm_group": "treatment"},
            {"sample_name": "sample3", "norm_group": "treatment"}
        ])
        return df
    
    design.to_dataframe = to_dataframe_with_groups
    return design


@pytest.fixture
def mock_design_ip():
    """Provide a mock DesignIP object."""
    return MockDesignIP(["exp1", "exp2"], ["chip1", "chip2"], ["input", "input"])


@pytest.mark.unit
class TestNormGroup:
    """Test cases for the NormGroup class."""

    def test_norm_group_creation(self):
        """Test basic NormGroup creation."""
        group = NormGroup(
            group="treatment",
            samples=["sample1", "sample2"],
            reference_sample="sample1"
        )
        assert group.group == "treatment"
        assert group.samples == ["sample1", "sample2"]
        assert group.reference_sample == "sample1"

    def test_norm_group_defaults(self):
        """Test NormGroup creation with defaults."""
        group = NormGroup(samples=["sample1"])
        assert group.group == "all"
        assert group.reference_sample is None

    def test_norm_group_len(self):
        """Test NormGroup length."""
        group = NormGroup(samples=["sample1", "sample2", "sample3"])
        assert len(group) == 3

    def test_norm_group_contains(self):
        """Test NormGroup membership testing."""
        group = NormGroup(samples=["sample1", "sample2"])
        assert "sample1" in group
        assert "sample2" in group
        assert "sample3" not in group

    def test_norm_group_str(self):
        """Test NormGroup string representation."""
        group = NormGroup(
            group="test", 
            samples=["s1", "s2"], 
            reference_sample="s1"
        )
        expected = "NormGroup(group='test', samples=2, reference='s1')"
        assert str(group) == expected

    def test_from_design_basic(self, mock_design):
        """Test creating NormGroup from basic design."""
        group = NormGroup.from_design(mock_design)
        assert group.group == "all"
        assert set(group.samples) == {"sample1", "sample2", "sample3"}
        assert group.reference_sample == "sample1"  # First sample

    def test_from_design_with_subset(self, mock_design_with_groups):
        """Test creating NormGroup with subset values."""
        group = NormGroup.from_design(
            mock_design_with_groups,
            subset_column="norm_group",
            subset_value=["treatment"]
        )
        assert group.group == "treatment"
        assert set(group.samples) == {"sample2", "sample3"}

    def test_from_design_ip_no_controls(self, mock_design_ip):
        """Test creating NormGroup from DesignIP without controls."""
        group = NormGroup.from_design(
            mock_design_ip,
            include_controls=False
        )
        # Should include IP samples with naming convention
        assert len(group.samples) > 0

    def test_from_design_ip_with_controls(self, mock_design_ip):
        """Test creating NormGroup from DesignIP with controls."""
        group = NormGroup.from_design(
            mock_design_ip,
            include_controls=True
        )
        # Should include both IP and control samples
        assert len(group.samples) > 0

    def test_from_design_empty_result_raises_error(self, mock_design):
        """Test that empty results raise ValueError."""
        with pytest.raises(ValueError, match="No samples found"):
            NormGroup.from_design(
                mock_design,
                subset_column="nonexistent",
                subset_value=["nonexistent"]
            )

    def test_prepare_dataframe_design(self, mock_design):
        """Test _prepare_dataframe with Design object."""
        df = NormGroup._prepare_dataframe(mock_design, include_controls=False)
        assert "sample_fullname" in df.index.names or df.index.name == "sample_fullname"
        assert len(df) == 3

    def test_prepare_dataframe_design_ip(self, mock_design_ip):
        """Test _prepare_dataframe with DesignIP object."""
        df = NormGroup._prepare_dataframe(mock_design_ip, include_controls=False)
        assert len(df) >= 2  # Should have IP samples


@pytest.mark.unit
class TestNormGroups:
    """Test cases for the NormGroups class."""

    def test_norm_groups_creation(self):
        """Test basic NormGroups creation."""
        group1 = NormGroup(group="control", samples=["sample1"])
        group2 = NormGroup(group="treatment", samples=["sample2", "sample3"])
        groups = NormGroups(groups=[group1, group2])
        
        assert len(groups) == 2
        assert len(groups.groups) == 2

    def test_from_design_no_groups_column(self, mock_design):
        """Test creating NormGroups when no grouping column exists."""
        groups = NormGroups.from_design(mock_design)
        assert len(groups) == 1
        assert groups.groups[0].group == "all"
        assert set(groups.groups[0].samples) == {"sample1", "sample2", "sample3"}

    def test_from_design_with_groups_column(self, mock_design_with_groups):
        """Test creating NormGroups with grouping column."""
        groups = NormGroups.from_design(
            mock_design_with_groups,
            subset_column="norm_group"
        )
        assert len(groups) == 2  # control and treatment
        
        group_names = {group.group for group in groups.groups}
        assert group_names == {"control", "treatment"}

    def test_sample_to_group_mapping(self):
        """Test sample to group mapping property."""
        group1 = NormGroup(group="A", samples=["sample1", "sample2"])
        group2 = NormGroup(group="B", samples=["sample3"])
        groups = NormGroups(groups=[group1, group2])
        
        mapping = groups.sample_to_group_mapping
        assert mapping["sample1"] == "A"
        assert mapping["sample2"] == "A"
        assert mapping["sample3"] == "B"

    def test_group_to_samples_mapping(self):
        """Test group to samples mapping property."""
        group1 = NormGroup(group="A", samples=["sample1", "sample2"])
        group2 = NormGroup(group="B", samples=["sample3"])
        groups = NormGroups(groups=[group1, group2])
        
        mapping = groups.group_to_samples_mapping
        assert mapping["A"] == ["sample1", "sample2"]
        assert mapping["B"] == ["sample3"]

    def test_get_sample_group(self):
        """Test getting group for a sample."""
        group1 = NormGroup(group="control", samples=["sample1"])
        group2 = NormGroup(group="treatment", samples=["sample2"])
        groups = NormGroups(groups=[group1, group2])
        
        assert groups.get_sample_group("sample1") == "control"
        assert groups.get_sample_group("sample2") == "treatment"

    def test_get_sample_group_not_found(self):
        """Test getting group for nonexistent sample raises KeyError."""
        group = NormGroup(group="test", samples=["sample1"])
        groups = NormGroups(groups=[group])
        
        with pytest.raises(KeyError, match="Sample 'nonexistent' not found"):
            groups.get_sample_group("nonexistent")

    def test_get_samples_in_group(self):
        """Test getting samples in a group."""
        group1 = NormGroup(group="A", samples=["sample1", "sample2"])
        group2 = NormGroup(group="B", samples=["sample3"])
        groups = NormGroups(groups=[group1, group2])
        
        assert groups.get_samples_in_group("A") == ["sample1", "sample2"]
        assert groups.get_samples_in_group("B") == ["sample3"]

    def test_get_samples_in_group_not_found(self):
        """Test getting samples for nonexistent group raises KeyError."""
        group = NormGroup(group="test", samples=["sample1"])
        groups = NormGroups(groups=[group])
        
        with pytest.raises(KeyError, match="Group 'nonexistent' not found"):
            groups.get_samples_in_group("nonexistent")

    def test_get_group_by_name(self):
        """Test getting NormGroup object by name."""
        group1 = NormGroup(group="control", samples=["sample1"])
        group2 = NormGroup(group="treatment", samples=["sample2"])
        groups = NormGroups(groups=[group1, group2])
        
        retrieved = groups.get_group_by_name("control")
        assert retrieved == group1
        assert retrieved.samples == ["sample1"]

    def test_get_group_by_name_not_found(self):
        """Test getting nonexistent group by name raises KeyError."""
        group = NormGroup(group="test", samples=["sample1"])
        groups = NormGroups(groups=[group])
        
        with pytest.raises(KeyError, match="Group 'nonexistent' not found"):
            groups.get_group_by_name("nonexistent")

    def test_groups_iteration(self):
        """Test iterating over groups."""
        group1 = NormGroup(group="A", samples=["sample1"])
        group2 = NormGroup(group="B", samples=["sample2"])
        groups = NormGroups(groups=[group1, group2])
        
        iterated = list(groups)
        assert len(iterated) == 2
        assert group1 in iterated
        assert group2 in iterated

    def test_groups_contains(self):
        """Test group membership testing."""
        group1 = NormGroup(group="control", samples=["sample1"])
        group2 = NormGroup(group="treatment", samples=["sample2"])
        groups = NormGroups(groups=[group1, group2])
        
        assert "control" in groups
        assert "treatment" in groups
        assert "nonexistent" not in groups

    def test_groups_str(self):
        """Test NormGroups string representation."""
        group1 = NormGroup(group="A", samples=["s1", "s2"])
        group2 = NormGroup(group="B", samples=["s3"])
        groups = NormGroups(groups=[group1, group2])
        
        result = str(groups)
        assert "NormGroups(2 groups:" in result
        assert "'A': 2 samples" in result
        assert "'B': 1 samples" in result

    # Test backward compatibility properties
    def test_backward_compatibility_sample_groups(self):
        """Test backward compatibility sample_groups property."""
        group1 = NormGroup(group="A", samples=["sample1"])
        groups = NormGroups(groups=[group1])
        
        # Should be alias for group_to_samples_mapping
        assert groups.sample_groups == groups.group_to_samples_mapping

    def test_backward_compatibility_group_samples(self):
        """Test backward compatibility group_samples property."""
        group1 = NormGroup(group="A", samples=["sample1"])
        groups = NormGroups(groups=[group1])
        
        # Should be alias for sample_to_group_mapping
        assert groups.group_samples == groups.sample_to_group_mapping

    def test_backward_compatibility_get_grouped_samples(self):
        """Test backward compatibility get_grouped_samples method."""
        group1 = NormGroup(group="A", samples=["sample1", "sample2"])
        groups = NormGroups(groups=[group1])
        
        # Should be alias for get_samples_in_group
        assert groups.get_grouped_samples("A") == groups.get_samples_in_group("A")


@pytest.mark.integration
class TestNormGroupsIntegration:
    """Integration tests for NormGroups with real design objects."""
    
    @pytest.fixture
    def sample_design_data(self, tmp_path):
        """Create sample design data for testing."""
        # This would require actual Design/DesignIP objects
        # For now, we'll skip this and focus on unit tests
        pytest.skip("Integration tests require full design objects")

    def test_integration_with_real_design(self, sample_design_data):
        """Test NormGroups with real Design objects."""
        # This would test the full integration
        pass


@pytest.mark.parametrize("group_name,samples,reference", [
    ("control", ["c1", "c2"], "c1"),
    ("treatment", ["t1", "t2", "t3"], "t1"),
    (1, ["sample1"], "sample1"),
    ("mixed", ["a", "b", "c", "d"], "b"),
])
def test_norm_group_parametrized_creation(group_name, samples, reference):
    """Parametrized test for NormGroup creation."""
    group = NormGroup(
        group=group_name,
        samples=samples,
        reference_sample=reference
    )
    assert group.group == group_name
    assert group.samples == samples
    assert group.reference_sample == reference
    assert len(group) == len(samples)


@pytest.mark.edge_case
class TestNormGroupsEdgeCases:
    """Test edge cases for NormGroup and NormGroups."""

    def test_empty_samples_list(self):
        """Test NormGroup with empty samples list."""
        group = NormGroup(samples=[])
        assert len(group) == 0
        assert "anything" not in group

    def test_single_sample_group(self):
        """Test NormGroup with single sample."""
        group = NormGroup(samples=["only_sample"])
        assert len(group) == 1
        assert "only_sample" in group

    def test_numeric_group_names(self):
        """Test NormGroup with numeric group names."""
        group = NormGroup(group=123, samples=["sample1"])
        assert group.group == 123
        assert isinstance(group.group, int)

    def test_empty_norm_groups(self):
        """Test NormGroups with no groups."""
        groups = NormGroups(groups=[])
        assert len(groups) == 0
        assert list(groups) == []
        assert groups.sample_to_group_mapping == {}
        assert groups.group_to_samples_mapping == {}

    def test_duplicate_samples_across_groups(self):
        """Test behavior with duplicate samples across groups."""
        group1 = NormGroup(group="A", samples=["sample1", "sample2"])
        group2 = NormGroup(group="B", samples=["sample1", "sample3"])  # sample1 duplicated
        groups = NormGroups(groups=[group1, group2])
        
        # The mapping should contain the last occurrence
        mapping = groups.sample_to_group_mapping
        # This documents current behavior - last group wins
        assert mapping["sample1"] == "B"
        assert mapping["sample2"] == "A"
        assert mapping["sample3"] == "B"
