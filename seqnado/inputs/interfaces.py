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
    """Detect file type by majority rule across a list of paths.

    Counts occurrences of known extensions (fastq, bam, bigwig). Returns the
    type with the highest count; if there's a tie or no known types, returns None.
    """
    import pandas as pd

    counts = {"fastq": 0, "bam": 0, "bigwig": 0}
    for p in paths:
        s = str(Path(p)).lower()
        if s.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
            counts["fastq"] += 1
        elif s.endswith(".bam"):
            counts["bam"] += 1
        elif s.endswith((".bigwig", ".bw")):
            counts["bigwig"] += 1

    counts = pd.Series(counts)
    counts = counts.sort_values(ascending=False)

    # Check that there's a clear winner
    # Check the that highest count is > 0 and strictly greater than the second
    if counts.iloc[0] > 0 and counts.iloc[0] > counts.iloc[1]:
        return counts.index[0]
    return None
