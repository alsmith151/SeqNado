"""BAM file collection classes.

These classes allow users to start analyses from aligned BAM files instead of
raw FASTQ data. A BAM file typically stores both R1 and R2 alignments so we
model a BAM as a single file per sample.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Iterable

from pydantic import BaseModel

from .core import Metadata, is_valid_path, BaseCollection
from seqnado import Assay


class BamFile(BaseModel):
    """Represents a single BAM file on disk."""

    path: Path

    def model_post_init(self, __context: dict[str, Any] | None) -> None:  # type: ignore[override]
        self.path = self.path.absolute()
        if not is_valid_path(self.path):  # pragma: no cover - simple guard
            raise FileNotFoundError(f"BAM file not found: {self.path}")

    @property
    def sample_id(self) -> str:
        """Infer sample id from filename stem (without .bam)."""
        return self.path.stem


class BamCollection(BaseCollection):
    """Collection of BAM files with optional per-sample metadata.

    Provides convenience constructors analogous to `SampleCollection` but
    without paired-end logic.
    """

    bam_files: list[BamFile]

    @property
    def primary_file_type(self) -> str:
        return "bam"

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------
    @property
    def sample_ids(self) -> list[str]:
        return [b.sample_id for b in self.bam_files]

    @property
    def sample_names(self) -> list[str]:  # maintain naming parity
        return self.sample_ids

    @property
    def bam_paths(self) -> list[Path]:
        return [b.path for b in self.bam_files]

    def get_file_paths(self, kind: str | None = None) -> list[Path]:
        if kind is None or kind == "bam":
            return self.bam_paths
        raise ValueError(f"Unsupported file kind '{kind}' for BamCollection")

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------
    @classmethod
    def from_files(
        cls,
        assay: Assay,
        files: Iterable[str | Path],
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **metadata_kwargs: Any,
    ) -> BamCollection:
        bam_files = [BamFile(path=Path(f)) for f in files]
        if not bam_files:
            raise ValueError("No BAM files provided")

        # Prepare metadata template if user passed simple kwargs
        metadata_param = BaseCollection._prepare_metadata_for_directory(
            metadata, **metadata_kwargs
        )

        metadata_list: list[Metadata] = [
            BaseCollection._build_metadata(b.sample_id, metadata_param, assay)
            for b in bam_files
        ]
        return cls(assay=assay, bam_files=bam_files, metadata=metadata_list)

    @classmethod
    def from_directory(
        cls,
        assay: Assay,
        directory: str | Path,
        glob_pattern: str = "*.bam",
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **metadata_kwargs: Any,
    ) -> BamCollection:
        dir_path = Path(directory)
        files = list(dir_path.rglob(glob_pattern))
        if not files:
            raise FileNotFoundError(f"No BAM files found in {dir_path}")
        return cls.from_files(
            assay=assay, files=files, metadata=metadata, **metadata_kwargs
        )

    @classmethod
    def from_dataframe(
        cls, assay: Assay, df: Any, validate_deseq2: bool = False, assay_for_validation: Assay | None = None, **kwargs: Any
    ) -> BamCollection:
        """Build a BamCollection from a DataFrame.
        
        Expects columns: sample_id, bam, plus any metadata fields.
        
        Args:
            assay: The assay type
            df: DataFrame with sample metadata
            validate_deseq2: If True, require deseq2 field to be non-null (for RNA assays)
            assay_for_validation: Assay type to check in validation context
            **kwargs: Additional kwargs
        """
        import pandas as pd
        
        bam_files: list[BamFile] = []
        metadata: list[Metadata] = []
        metadata_fields = set(Metadata.model_fields.keys())
        
        # Use provided assay_for_validation or fall back to assay
        validation_assay = assay_for_validation or assay

        for rec in df.to_dict(orient="records"):
            # Build BamFile
            bam_files.append(BamFile(path=Path(rec["bam"])))

            # Collect metadata with validation context
            meta_fields = {k: rec.get(k) for k in metadata_fields if k in rec}
            metadata.append(Metadata.model_validate(meta_fields, context={'validate_deseq2': validate_deseq2, 'assay': validation_assay}))

        return cls(assay=assay, bam_files=bam_files, metadata=metadata)

    # ------------------------------------------------------------------
    # Data export
    # ------------------------------------------------------------------
    def to_dataframe(self):  # return pd.DataFrame but keep optional import
        import pandas as pd

        rows: list[dict[str, Any]] = []
        for bam, md in zip(self.bam_files, self.metadata):
            row = {"sample_id": bam.sample_id, "bam": bam.path}
            row.update(md.model_dump(exclude_none=True))
            rows.append(row)
        return pd.DataFrame(rows).sort_values("sample_id")

    # ------------------------------------------------------------------
    # File operations
    # ------------------------------------------------------------------
    def symlink_bam_files(self, output_dir: str | Path) -> None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        for bam in self.bam_files:
            link = out / f"{bam.sample_id}.bam"
            if not link.exists():
                link.symlink_to(bam.path.resolve())
