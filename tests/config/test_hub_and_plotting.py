from pathlib import Path
import pytest

from seqnado import GenomicCoordinate
from seqnado.config import UCSCHubConfig


def test_ucsc_hub_validation_and_fields(tmp_path: Path):
    # Construct explicitly and check fields validate
    parent = tmp_path / "parent"
    parent.mkdir()
    hub = UCSCHubConfig(
        directory=str(parent / "hub"),
        name="seqnado_hub",
        genome="hg38",
        email="a@b.com",
        default_position=GenomicCoordinate(chromosome="chr1", start=100, end=200),
        subgroup_by=["method", "norm", "strand"],
        overlay_by=["samplename", "method", "norm"],
    )
    assert "strand" in hub.subgroup_by
    assert set(["samplename", "method", "norm"]).issubset(hub.overlay_by)

    # Name validation and directory parent existence
    with pytest.raises(ValueError):
        UCSCHubConfig(
            directory=str(tmp_path / "missing" / "hub"),
            name="ok",
            default_position=GenomicCoordinate.from_string("chr1:100-200"),
        )
    ok = UCSCHubConfig(directory=str(parent / "hub"), name="seqnado_hub")
    assert ok.name == "seqnado_hub"
    with pytest.raises(ValueError):
        UCSCHubConfig(directory=str(parent / "hub"), name="bad name")
