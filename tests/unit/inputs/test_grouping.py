"""Tests for seqnado.inputs.grouping module."""

import pandas as pd
import pytest

from seqnado.inputs.grouping import SampleGroup, SampleGroupings, SampleGroups


class TestSampleGroup:
    """Tests for SampleGroup class."""

    def test_create_sample_group(self):
        """Test creating a basic sample group."""
        group = SampleGroup(name="batch1", samples=["s1", "s2", "s3"])
        
        assert group.name == "batch1"
        assert len(group.samples) == 3
        assert group.reference_sample is None

    def test_sample_group_with_reference(self):
        """Test creating a sample group with a reference sample."""
        group = SampleGroup(
            name="batch1",
            samples=["s1", "s2", "s3"],
            reference_sample="s1"
        )
        
        assert group.reference_sample == "s1"

    def test_sample_group_length(self):
        """Test __len__ method."""
        group = SampleGroup(name="batch1", samples=["s1", "s2", "s3"])
        assert len(group) == 3

    def test_sample_group_contains(self):
        """Test __contains__ method."""
        group = SampleGroup(name="batch1", samples=["s1", "s2", "s3"])
        assert "s1" in group
        assert "s4" not in group

    def test_sample_group_str(self):
        """Test __str__ method."""
        group = SampleGroup(
            name="batch1",
            samples=["s1", "s2"],
            reference_sample="s1"
        )
        result = str(group)
        assert "batch1" in result
        assert "2 samples" in result
        assert "ref=s1" in result


class TestSampleGroups:
    """Tests for SampleGroups class."""

    def test_create_sample_groups(self):
        """Test creating SampleGroups with multiple groups."""
        g1 = SampleGroup(name="batch1", samples=["s1", "s2"])
        g2 = SampleGroup(name="batch2", samples=["s3", "s4"])
        groups = SampleGroups(groups=[g1, g2])
        
        assert len(groups) == 2

    def test_from_dataframe_basic(self):
        """Test creating SampleGroups from a DataFrame."""
        df = pd.DataFrame({
            "sample": ["s1", "s2", "s3", "s4"],
            "scaling_group": ["batch1", "batch1", "batch2", "batch2"],
        })
        df = df.set_index("sample")
        
        groups = SampleGroups.from_dataframe(df, subset_column="scaling_group")
        
        assert len(groups) == 2
        assert groups.get_samples("batch1") == ["s1", "s2"]
        assert groups.get_samples("batch2") == ["s3", "s4"]

    def test_from_dataframe_with_reference_sample(self):
        """Test creating SampleGroups with a specific reference sample."""
        df = pd.DataFrame({
            "sample": ["s1", "s2", "s3"],
            "scaling_group": ["batch1", "batch1", "batch1"],
        })
        df = df.set_index("sample")
        
        groups = SampleGroups.from_dataframe(
            df,
            subset_column="scaling_group",
            reference_sample="s2"
        )
        
        group = groups.get_group("batch1")
        assert group.reference_sample == "s2"

    def test_from_dataframe_missing_column_returns_empty(self):
        """Test that missing column returns empty SampleGroups."""
        df = pd.DataFrame({
            "sample": ["s1", "s2"],
            "other_column": ["a", "b"],
        })
        df = df.set_index("sample")
        
        groups = SampleGroups.from_dataframe(df, subset_column="missing_column")
        
        assert groups.empty
        assert len(groups) == 0

    def test_sample_to_group_mapping(self):
        """Test sample_to_group method."""
        g1 = SampleGroup(name="batch1", samples=["s1", "s2"])
        g2 = SampleGroup(name="batch2", samples=["s3", "s4"])
        groups = SampleGroups(groups=[g1, g2])
        
        mapping = groups.sample_to_group
        
        assert mapping["s1"] == "batch1"
        assert mapping["s2"] == "batch1"
        assert mapping["s3"] == "batch2"
        assert mapping["s4"] == "batch2"

    def test_group_to_samples_mapping(self):
        """Test group_to_samples method."""
        g1 = SampleGroup(name="batch1", samples=["s1", "s2"])
        g2 = SampleGroup(name="batch2", samples=["s3", "s4"])
        groups = SampleGroups(groups=[g1, g2])
        
        mapping = groups.group_to_samples
        
        assert mapping["batch1"] == ["s1", "s2"]
        assert mapping["batch2"] == ["s3", "s4"]

    def test_get_group_existing(self):
        """Test getting an existing group by name."""
        g1 = SampleGroup(name="batch1", samples=["s1", "s2"])
        groups = SampleGroups(groups=[g1])
        
        group = groups.get_group("batch1")
        assert group.name == "batch1"

    def test_get_group_missing_raises_error(self):
        """Test that getting a non-existent group raises KeyError."""
        groups = SampleGroups(groups=[])
        
        with pytest.raises(KeyError, match="Group 'missing' not found"):
            groups.get_group("missing")

    def test_get_samples_for_group(self):
        """Test getting samples for a specific group."""
        g1 = SampleGroup(name="batch1", samples=["s1", "s2", "s3"])
        groups = SampleGroups(groups=[g1])
        
        samples = groups.get_samples("batch1")
        assert samples == ["s1", "s2", "s3"]

    def test_empty_property_false(self):
        """Test empty property when groups exist."""
        g1 = SampleGroup(name="batch1", samples=["s1"])
        groups = SampleGroups(groups=[g1])
        
        assert not groups.empty

    def test_empty_property_true(self):
        """Test empty property when no groups exist."""
        groups = SampleGroups(groups=[])
        
        assert groups.empty

    def test_str_representation(self):
        """Test __str__ method."""
        g1 = SampleGroup(name="batch1", samples=["s1"])
        g2 = SampleGroup(name="batch2", samples=["s2"])
        groups = SampleGroups(groups=[g1, g2])
        
        result = str(groups)
        assert "2 groups" in result
        assert "batch1" in result
        assert "batch2" in result


class TestSampleGroupings:
    """Tests for SampleGroupings class."""

    def test_create_empty_sample_groupings(self):
        """Test creating empty SampleGroupings."""
        groupings = SampleGroupings()
        
        assert groupings.empty
        assert len(groupings.groupings) == 0

    def test_add_grouping(self):
        """Test adding a grouping to SampleGroupings."""
        groupings = SampleGroupings()
        g1 = SampleGroup(name="batch1", samples=["s1", "s2"])
        groups = SampleGroups(groups=[g1])
        
        groupings.add_grouping("scaling", groups)
        
        assert "scaling" in groupings
        assert not groupings.empty

    def test_get_grouping_existing(self):
        """Test getting an existing grouping."""
        groupings = SampleGroupings()
        g1 = SampleGroup(name="batch1", samples=["s1", "s2"])
        groups = SampleGroups(groups=[g1])
        groupings.add_grouping("scaling", groups)
        
        result = groupings.get_grouping("scaling")
        assert result == groups

    def test_get_grouping_missing_raises_error(self):
        """Test that getting a non-existent grouping raises KeyError."""
        groupings = SampleGroupings()
        
        with pytest.raises(KeyError, match="Grouping 'missing' not found"):
            groupings.get_grouping("missing")

    def test_contains_operator(self):
        """Test __contains__ operator."""
        groupings = SampleGroupings()
        g1 = SampleGroup(name="batch1", samples=["s1"])
        groups = SampleGroups(groups=[g1])
        groupings.add_grouping("scaling", groups)
        
        assert "scaling" in groupings
        assert "normalization" not in groupings

    def test_str_representation(self):
        """Test __str__ method."""
        groupings = SampleGroupings()
        g1 = SampleGroup(name="batch1", samples=["s1"])
        groups = SampleGroups(groups=[g1])
        groupings.add_grouping("scaling", groups)
        groupings.add_grouping("normalization", groups)
        
        result = str(groupings)
        assert "scaling" in result
        assert "normalization" in result

    def test_all_samples(self):
        """Test all_samples method returns unique samples."""
        groupings = SampleGroupings()
        
        g1 = SampleGroup(name="batch1", samples=["s1", "s2"])
        g2 = SampleGroup(name="batch2", samples=["s2", "s3"])
        groups1 = SampleGroups(groups=[g1])
        groups2 = SampleGroups(groups=[g2])
        
        groupings.add_grouping("scaling", groups1)
        groupings.add_grouping("normalization", groups2)
        
        all_samples = groupings.all_samples()
        assert set(all_samples) == {"s1", "s2", "s3"}

    def test_empty_property_false(self):
        """Test empty property when groupings exist."""
        groupings = SampleGroupings()
        g1 = SampleGroup(name="batch1", samples=["s1"])
        groups = SampleGroups(groups=[g1])
        groupings.add_grouping("scaling", groups)
        
        assert not groupings.empty

    def test_empty_property_true(self):
        """Test empty property when no groupings exist."""
        groupings = SampleGroupings()
        
        assert groupings.empty

    def test_multiple_groupings_independent(self):
        """Test that multiple groupings can be added independently."""
        groupings = SampleGroupings()
        
        g1 = SampleGroup(name="batch1", samples=["s1", "s2"])
        g2 = SampleGroup(name="condition1", samples=["s1", "s3"])
        
        groups_scaling = SampleGroups(groups=[g1])
        groups_condition = SampleGroups(groups=[g2])
        
        groupings.add_grouping("scaling", groups_scaling)
        groupings.add_grouping("condition", groups_condition)
        
        assert groupings.get_grouping("scaling").get_samples("batch1") == ["s1", "s2"]
        assert groupings.get_grouping("condition").get_samples("condition1") == ["s1", "s3"]
