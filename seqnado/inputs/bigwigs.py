"""BigWig file collection classes.

Supports grouping +/- strand tracks for RNA assays. For non-RNA assays one
bigWig per sample is expected.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Iterable

from pydantic import BaseModel

from .core import Metadata, is_valid_path, BaseCollection
from seqnado import Assay


class BigWigFile(BaseModel):
	path: Path

	def model_post_init(self, __context: dict[str, Any] | None) -> None:  # type: ignore[override]
		self.path = self.path.absolute()
		if not is_valid_path(self.path):  # pragma: no cover
			raise FileNotFoundError(f"bigWig file not found: {self.path}")

	@property
	def stem(self) -> str:
		return self.path.stem

	def infer_sample_id(self, assay: Assay) -> str:
		stem = self.stem
		if assay == Assay.RNA and (stem.endswith("_plus") or stem.endswith("_minus")):
			return stem.rsplit("_", 1)[0]
		return stem

	def is_strand_specific(self, assay: Assay) -> bool:
		stem = self.stem
		return assay == Assay.RNA and (stem.endswith("_plus") or stem.endswith("_minus"))

	def strand(self, assay: Assay) -> str | None:
		if not self.is_strand_specific(assay):
			return None
		return self.stem.rsplit("_", 1)[1]


class BigWigCollection(BaseModel):
	assay: Assay
	bigwig_files: list[BigWigFile]
	metadata: list[Metadata]

	@property
	def primary_file_type(self) -> str:
		return "bigwig"

	# ---------------------------------------------------------------
	# Properties
	# ---------------------------------------------------------------
	@property
	def sample_ids(self) -> list[str]:
		return [bw.infer_sample_id(self.assay) for bw in self.bigwig_files]

	@property
	def sample_names(self) -> list[str]:  # parity with SampleCollection
		return sorted(set(self.sample_ids))

	@property
	def tracks(self) -> dict[str, list[Path]]:
		mapping: dict[str, list[Path]] = {}
		for bw in self.bigwig_files:
			sid = bw.infer_sample_id(self.assay)
			mapping.setdefault(sid, []).append(bw.path)
		# ensure deterministic ordering: plus before minus if RNA
		for sid, paths in mapping.items():
			mapping[sid] = sorted(
				paths, key=lambda p: ("minus" in p.stem, p.name)
			)
		return mapping

	def get_file_paths(self, kind: str | None = None) -> list[Path]:
		if kind is None or kind == "bigwig":
			return [bw.path for bw in self.bigwig_files]
		raise ValueError(f"Unsupported file kind '{kind}' for BigWigCollection")

	
	# ---------------------------------------------------------------
	# Constructors
	# ---------------------------------------------------------------
	@classmethod
	def from_files(
		cls,
		assay: Assay,
		files: Iterable[str | Path],
		metadata: Callable[[str], Metadata] | Metadata | None = None,
		**metadata_kwargs: Any,
	) -> BigWigCollection:
		bigwig_files = [BigWigFile(path=Path(f)) for f in files]
		if not bigwig_files:
			raise ValueError("No bigWig files provided")

		metadata_param = BaseCollection._prepare_metadata_for_directory(
			metadata, **metadata_kwargs
		)

		metadata_list: list[Metadata] = []
		for bw in bigwig_files:
			sid = bw.infer_sample_id(assay)
			metadata_list.append(
				BaseCollection._build_metadata(sid, metadata_param, assay)
			)

		return cls(assay=assay, bigwig_files=bigwig_files, metadata=metadata_list)

	@classmethod
	def from_directory(
		cls,
		assay: Assay,
		directory: str | Path,
		glob_patterns: tuple[str, ...] = ("*.bigWig", "*.bw"),
		metadata: Callable[[str], Metadata] | Metadata | None = None,
		**metadata_kwargs: Any,
	) -> BigWigCollection:
		dir_path = Path(directory)
		files: list[Path] = []
		for pattern in glob_patterns:
			files.extend(dir_path.rglob(pattern))
		if not files:
			raise FileNotFoundError(f"No bigWig files found in {dir_path}")
		return cls.from_files(
			assay=assay, files=files, metadata=metadata, **metadata_kwargs
		)

	# ---------------------------------------------------------------
	# Export
	# ---------------------------------------------------------------
	def to_dataframe(self):  # returns pd.DataFrame
		import pandas as pd

		# group tracks per sample
		rows: list[dict[str, Any]] = []
		for sid, paths in self.tracks.items():
			row: dict[str, Any] = {"sample_id": sid}
			if self.assay == Assay.RNA:
				plus = next((p for p in paths if p.stem.endswith("_plus")), None)
				minus = next((p for p in paths if p.stem.endswith("_minus")), None)
				row.update({"bigwig_plus": plus, "bigwig_minus": minus})
			else:
				row.update({"bigwig": paths[0]})

			# pick metadata corresponding to first path sample
			md_index = self.sample_ids.index(sid)
			row.update(self.metadata[md_index].model_dump(exclude_none=True))
			rows.append(row)
		return pd.DataFrame(rows).sort_values("sample_id")

	# ---------------------------------------------------------------
	# File ops
	# ---------------------------------------------------------------
	def symlink_bigwig_files(self, output_dir: str | Path) -> None:
		out = Path(output_dir)
		out.mkdir(parents=True, exist_ok=True)
		for bw in self.bigwig_files:
			link = out / bw.path.name
			if not link.exists():
				link.symlink_to(bw.path.resolve())

