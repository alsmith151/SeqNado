import importlib
from enum import Enum
import pytest
from pydantic import BaseModel
from seqnado.inputs.interfaces import detect_file_type

# language: python
# filepath: /workspaces/SeqNado/seqnado/inputs/test___init__.py


# Relative imports from the package under test
from seqnado.inputs import (
    Assay,
    Metadata,
    FastqCollection,
    FastqCollectionForIP,
    SampleGroups,
    SampleGroupings,
    BamFile,
    BamCollection,
    BigWigFile,
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
    assert expected.issubset(set(inputs_all))


def test_module_attributes_exist_for___all__():
    pkg = importlib.import_module("seqnado.inputs")
    for name in inputs_all:
        assert hasattr(pkg, name), f"Missing attribute {name} on seqnado.inputs"


def test_assay_is_enum_and_has_members():
    assert issubclass(Assay, Enum)
    assert len(list(Assay)) >= 1


def test_metadata_is_pydantic_model_and_instantiable():
    assert issubclass(Metadata, BaseModel)
    m = Metadata()  # should use defaults
    assert isinstance(m, Metadata)


def test_collectionlike_protocol_accepts_minimal_impl_and_rejects_invalid():
    # Minimal object implementing the protocol
    class MinimalCollection:
        def __init__(self):
            # pick first available Assay member if possible, else fallback to attribute presence
            self.assay = next(iter(Assay)) if len(list(Assay)) else "dummy"

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
    ],
)
def test_detect_file_type_heuristics(paths, expected):
    assert detect_file_type(paths) == expected