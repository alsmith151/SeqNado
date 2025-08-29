"""Interfaces / Protocols for input collections.

Provides a common, minimal surface so the pipeline can operate on any
collection type (FASTQ, IP FASTQ, BAM, bigWig) without hard unions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol, Iterable, runtime_checkable

from seqnado import Assay


@runtime_checkable
class CollectionLike(Protocol):
    """A minimal protocol implemented by all input collection types."""

    assay: Assay

    @property
    def sample_names(self) -> list[str]:
        ...  # pragma: no cover - interface only

    def to_dataframe(self):  # returns a pandas DataFrame; typing avoided to keep optional dep
        ...  # pragma: no cover

    @property
    def primary_file_type(self) -> str:
        """Primary file type backing this collection (e.g. fastq, bam, bigwig)."""
        ...  # pragma: no cover

    def get_file_paths(self, kind: str | None = None) -> list[Path]:
        """Return list of underlying file paths.

        kind defaults to the collection's primary file type. Implementations may
        accept additional granular kinds (e.g. fastq_r1, fastq_r2).
        """
        ...  # pragma: no cover


def ensure_collection(obj) -> CollectionLike:
    """Validate an arbitrary object implements CollectionLike."""
    if not isinstance(obj, CollectionLike):  # runtime_checkable enables this
        raise TypeError(
            "Object does not implement required collection interface: CollectionLike"
        )
    return obj


def detect_file_type(paths: Iterable[str | Path]) -> str | None:
    """Simple heuristic to detect file type from a list of paths."""
    exts = {str(Path(p)).lower() for p in paths}
    if any(e.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')) for e in exts):
        return 'fastq'
    if any(e.endswith('.bam') for e in exts):
        return 'bam'
    if any(e.endswith(('.bigwig', '.bw')) for e in exts):
        return 'bigwig'
    return None
