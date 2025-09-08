import importlib

from seqnado.inputs import __all__ as inputs_all


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
