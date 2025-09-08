from enum import Enum

import pytest
from pydantic import BaseModel, ValidationError

from seqnado.inputs import Assay, CollectionLike, Metadata, ensure_collection


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


def test_metadata_validator_rejects_none_values():
    with pytest.raises(ValidationError):
        Metadata(consensus_group=None)
