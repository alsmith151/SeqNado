import pandas as pd
from seqnado.inputs.grouping import SampleGroup, SampleGroups, SampleGroupings


def test_samplegroups_from_dataframe():
    df = pd.DataFrame({
        "scaling_group": ["g1", "g1", "g2"],
    }, index=["s1", "s2", "s3"])

    groups = SampleGroups.from_dataframe(df, subset_column="scaling_group")
    assert len(groups) == 2
    names = {g.name for g in groups.groups}
    assert names == {"g1", "g2"}
    assert set(groups.get_samples("g1")) == {"s1", "s2"}


def test_samplegroupings_container():
    g1 = SampleGroups(groups=[SampleGroup(name="g", samples=["a", "b"])])
    sg = SampleGroupings()
    sg.add_grouping("consensus", g1)
    assert "consensus" in sg
    assert sg.get_grouping("consensus").groups[0].name == "g"