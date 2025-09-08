import importlib
from pathlib import Path
from enum import Enum
import pytest
from pydantic import ValidationError
from pydantic import BaseModel
from seqnado.inputs.interfaces import detect_file_type

# language: python
# filepath: /workspaces/SeqNado/tests/test_inputs_submodule.py


# Import the function under test via a module-relative path

# Import package public API
from seqnado.inputs import (
    Assay,
    Metadata,
    BamCollection,
    BigWigCollection,
    CollectionLike,
    ensure_collection,
    __all__ as inputs_all,
)


def test_public_api___all__contains_expected_exports():
    expected = {
        "Assay",
        "Metadata",
        "FastqCollection",
        "FastqCollectionForIP",
        "SampleGroups",
        "SampleGroupings",
        "BamFile",
        "BamCollection",
        "BigWigFile",
        "BigWigCollection",
        "CollectionLike",
        "ensure_collection",
    }
    assert isinstance(inputs_all, list)
    assert expected.issubset(set(inputs_all))


def test_module_attributes_exist_for___all__():
    pkg = importlib.import_module("seqnado.inputs")
    for name in inputs_all:
        assert hasattr(pkg, name), f"Missing attribute {name} on seqnado.inputs"


def test___all___has_unique_entries_and_no_privates():
    assert len(inputs_all) == len(set(inputs_all)), "__all__ contains duplicates"
    assert all(not n.startswith("_") for n in inputs_all), "__all__ exposes private names"


def test_assay_is_enum_and_has_members():
    assert issubclass(Assay, Enum)
    assert len(list(Assay)) >= 1


def test_metadata_is_pydantic_model_and_instantiable():
    assert issubclass(Metadata, BaseModel)
    m = Metadata()  # should construct with defaults
    assert isinstance(m, Metadata)


def test_collectionlike_protocol_accepts_minimal_impl_and_rejects_invalid():
    class MinimalCollection:
        def __init__(self):
            self.assay = next(iter(Assay)) if len(list(Assay)) else None

        @property
        def sample_names(self):
            return []

        def to_dataframe(self):
            return None

        @property
        def primary_file_type(self):
            return "fastq"

        def get_file_paths(self, kind=None):
            return []

    ok = MinimalCollection()
    assert isinstance(ok, CollectionLike)
    assert ensure_collection(ok) is ok

    class NotACollection:
        pass

    with pytest.raises(TypeError):
        ensure_collection(NotACollection())


@pytest.mark.parametrize(
    "paths,expected",
    [
        (["a.fastq", "b.fastq.gz", "c.FQ"], "fastq"),
        (["x.bam"], "bam"),
        (["x.bigWig", "y.bw"], "bigwig"),
        (["unknown.txt", "notes.md"], None),
        # Majority rule
        (["a.fastq", "b.bam", "c.fastq.gz"], "fastq"),
        (["a.bam", "b.bam", "c.bw"], "bam"),
        # Tie -> None
        (["a.fastq", "b.bam"], None),
        (["a.bw", "b.bigWig"], "bigwig"),
    ],
)
def test_detect_file_type_majority(paths, expected):
    assert detect_file_type(paths) == expected


# ---------------------------------
# Advanced example-based tests
# ---------------------------------

def _touch(p: Path) -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("")
    return p


def test_bamcollection_from_files_and_dataframe(tmp_path: Path):
    f1 = _touch(tmp_path / "a.bam")
    f2 = _touch(tmp_path / "b.bam")

    def md_builder(sample: str):
        from seqnado.inputs import Metadata

        return Metadata(consensus_group="g1")

    coll = BamCollection.from_files(
        assay=Assay.ATAC, files=[f1, f2], metadata=md_builder
    )
    assert coll.primary_file_type == "bam"
    assert sorted(coll.sample_ids) == ["a", "b"]
    assert set(coll.get_file_paths()) == {f1, f2}
    df = coll.to_dataframe()
    # Expect rows for both samples and metadata columns present
    assert set(df["sample_id"]) == {"a", "b"}
    assert "consensus_group" in df.columns and "assay" in df.columns


def test_bamcollection_get_file_paths_kinds_and_symlinks(tmp_path: Path):
    f1 = _touch(tmp_path / "s1.bam")
    coll = BamCollection.from_files(assay=Assay.ATAC, files=[f1])
    # kind filter
    assert coll.get_file_paths("bam") == [f1]
    with pytest.raises(ValueError):
        coll.get_file_paths("fastq")
    # symlink
    out = tmp_path / "links"
    coll.symlink_bam_files(out)
    assert (out / "s1.bam").exists()


def test_bigwigcollection_rna_strand_logic_and_dataframe(tmp_path: Path):
    s1_plus = _touch(tmp_path / "s1_plus.bw")
    s1_minus = _touch(tmp_path / "s1_minus.bw")
    s2_plus = _touch(tmp_path / "s2_plus.bigWig")
    s2_minus = _touch(tmp_path / "s2_minus.bigWig")

    coll = BigWigCollection.from_files(
        assay=Assay.RNA,
        files=[s1_minus, s2_plus, s1_plus, s2_minus],
    )
    # sample ids derive from stem without strand suffix
    assert sorted(coll.sample_ids) == ["s1", "s1", "s2", "s2"]
    # tracks map to lists ordered plus before minus
    tracks = coll.tracks
    assert [p.stem for p in tracks["s1"]] == ["s1_plus", "s1_minus"]
    assert [p.stem for p in tracks["s2"]] == ["s2_plus", "s2_minus"]
    # dataframe has plus/minus columns
    df = coll.to_dataframe()
    assert set(df["sample_id"]) == {"s1", "s2"}
    assert {"bigwig_plus", "bigwig_minus"}.issubset(df.columns)


def test_bigwigcollection_non_rna_single_track_and_symlinks(tmp_path: Path):
    s1 = _touch(tmp_path / "s1.bw")
    s2 = _touch(tmp_path / "s2.bigWig")
    coll = BigWigCollection.from_files(assay=Assay.ATAC, files=[s1, s2])
    assert coll.primary_file_type == "bigwig"
    assert sorted(coll.sample_names) == ["s1", "s2"]
    df = coll.to_dataframe()
    assert "bigwig" in df.columns and set(df["sample_id"]) == {"s1", "s2"}
    out = tmp_path / "bw_links"
    coll.symlink_bigwig_files(out)
    assert (out / s1.name).exists() and (out / s2.name).exists()


def test_ensure_collection_accepts_bam_and_bigwig(tmp_path: Path):
    b = BamCollection.from_files(assay=Assay.ATAC, files=[_touch(tmp_path / "a.bam")])
    bw = BigWigCollection.from_files(assay=Assay.RNA, files=[_touch(tmp_path / "s_plus.bw")])
    assert ensure_collection(b) is b
    assert ensure_collection(bw) is bw


def test_metadata_validator_rejects_none_values():
    from seqnado.inputs import Metadata

    with pytest.raises(ValidationError):
        Metadata(consensus_group=None)