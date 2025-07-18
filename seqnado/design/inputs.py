# refactored_design.py
from __future__ import annotations
import pathlib
from collections import defaultdict
from itertools import chain
from typing import Any, Callable, Iterable, Literal
import pandas as pd
from pydantic import BaseModel, field_validator
from pandera.typing import DataFrame, Series
from .fastq import FastqFile, FastqSet, FastqSetIP, FastqFileIP
from .core import Metadata, Assay
from .validation import DesignDataFrame
from .experiment import IPExperiment


class Design(BaseModel):
    """
    Represents a collection of sequencing samples (FASTQ files) grouped into named sets,
    with optional per-sample metadata.

    Attributes:
        fastq_sets: List of FastqSet objects (paired or single-end samples).
        metadata:  List of Metadata objects corresponding one-to-one with fastq_sets.
    """
    assay: Assay
    fastq_sets: list[FastqSet]
    metadata: list[Metadata]

    @field_validator("assay")
    @classmethod
    def validate_non_ip_assay(cls, v: Assay) -> Assay:
        """Ensure the assay doesn't require IP (immunoprecipitation)."""
        if v in Assay.ip_assays():
            raise ValueError(f"Assay '{v.value}' requires IP and should use DesignIP instead")
        return v

    @property
    def sample_names(self) -> list[str]:
        """
        Returns all sample names in the design.
        """
        return [fs.name for fs in self.fastq_sets]

    @property
    def fastq_paths(self) -> list[pathlib.Path]:
        """
        Flattens all R1/R2 file paths into a single list.
        """
        return [
            path
            for fs in self.fastq_sets
            for path in (fs.r1.path, *(fs.r2.path if fs.r2 else []))
        ]

    def query(self, sample_name: str) -> FastqSet:
        """
        Retrieve the FastqSet by its sample name.

        Raises:
            ValueError: If sample_name not found.
        """
        try:
            return next(fs for fs in self.fastq_sets if fs.name == sample_name)
        except StopIteration:
            raise ValueError(f"Sample '{sample_name}' not found in Design")

    @classmethod
    def from_fastq_files(
        cls,
        assay: Assay,
        files: Iterable[str | pathlib.Path],
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **fastqset_kwargs: Any,
    ) -> Design:
        """
        Build a Design by scanning a list of FASTQ paths:

        1. Convert raw paths to FastqFile.
        2. Group by `sample_base` and sort by read_number.
        3. Create FastqSet (single- or paired-end) for each sample.
        4. Generate Metadata via `metadata(sample_name)`, or default.

        Args:
            files: Iterable of file paths (strings or pathlib.Path).
            metadata:
                - Callable(sample_name) → Metadata to customize per-sample metadata.
                - Single Metadata instance applied to all.
                - None → defaults to Metadata(scale_group="all").
            fastqset_kwargs: Extra fields forwarded to FastqSet constructor.
        """
        # Convert and sort
        fq_files = [FastqFile(path=pathlib.Path(f)) for f in files]
        fq_files.sort(key=lambda x: (x.sample_base, x.read_number))

        # Group by sample_stem
        groups: dict[str, list[FastqFile]] = defaultdict(list)
        for fq in fq_files:
            groups[fq.sample_base].append(fq)

        _fastq_sets: list[FastqSet] = []
        _metadata: list[Metadata] = []
        for sample, fqs in groups.items():
            # Build FastqSet
            if len(fqs) == 1:
                fs = FastqSet(name=sample, r1=fqs[0], **fastqset_kwargs)
            elif len(fqs) == 2:
                fs = FastqSet(name=sample, r1=fqs[0], r2=fqs[1], **fastqset_kwargs)
            else:
                raise ValueError(
                    f"Unexpected number of FASTQ files for '{sample}': {len(fqs)}"
                )
            _fastq_sets.append(fs)

            # Build Metadata
            if callable(metadata):
                meta = metadata(sample)
            elif isinstance(metadata, Metadata):
                meta = metadata
            else:
                meta = Metadata(scale_group="all")
            _metadata.append(meta)

        return cls(assay=assay, fastq_sets=_fastq_sets, metadata=_metadata)

    @classmethod
    def from_directory(
        cls,
        assay: Assay,
        directory: str | pathlib.Path,
        glob_patterns: Iterable[str] = ("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"),
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **kwargs: Any,
    ) -> Design:
        """
        Recursively scan a directory for FASTQ files and build a Design.

        Args:
            directory: Root path to search.
            glob_patterns: Filename patterns to include.
            metadata: Callable(sample_name) → Metadata or single Metadata instance.
            **kwargs: Extra fields converted directly to a shared Metadata.
        """
        dir_path = pathlib.Path(directory)
        files = list(chain.from_iterable(dir_path.glob(p) for p in glob_patterns))
        if not files:
            raise FileNotFoundError(f"No FASTQ files found in {dir_path}")

        if not callable(metadata) and not isinstance(metadata, Metadata):
            # If metadata is not callable, assume it's a single Metadata instance and add the kwargs to it
            metadata = Metadata(**{"scale_group": "all", **kwargs})

        return cls.from_fastq_files(assay=assay, files=files, metadata=metadata)

    def to_dataframe(self) -> pd.DataFrame:
        """
        Export the design to a pandas DataFrame, validated by DataFrameDesign.

        Columns: sample_name, r1, r2, plus all metadata fields.
        """
        import pandas as pd

        rows: list[dict[str, Any]] = []
        for fs, md in zip(self.fastq_sets, self.metadata):
            row: dict[str, Any] = {
                "sample_name": fs.name,
                "r1": fs.r1.path,
                "r2": fs.r2.path if fs.r2 else None,
            }
            metadata_dict = md.model_dump(exclude_none=True)
            # Convert Assay enum to string value for schema validation
            if 'assay' in metadata_dict and hasattr(metadata_dict['assay'], 'value'):
                metadata_dict['assay'] = metadata_dict['assay'].value
            row.update(metadata_dict)
            rows.append(row)

        df = pd.DataFrame(rows).sort_values("sample_name")
        return DataFrame[DesignDataFrame](df)

    @classmethod
    def from_dataframe(cls, assay: Assay, df: pd.DataFrame, **fastqset_kwargs: Any) -> Design:
        """
        Build a Design from a DataFrame, validated by DataFrameDesign.

        Expects columns: sample_name, r1, r2, plus any metadata fields.
        """
        df = DesignDataFrame.validate(df)
        fastq_sets: list[FastqSet] = []
        metadata: list[Metadata] = []
        metadata_fields = set(Metadata.model_fields.keys())

        for rec in df.to_dict(orient="records"):
            # Build FastqSet
            r2_path = rec.get("r2")
            fs = FastqSet(
                name=rec["sample_name"],
                r1=FastqFile(path=rec["r1"]),
                r2=FastqFile(path=r2_path) if r2_path else None,
                **fastqset_kwargs,
            )
            fastq_sets.append(fs)

            # Collect metadata
            meta_fields = {k: rec.get(k) for k in metadata_fields if k in rec}
            metadata.append(Metadata(**meta_fields))

        return cls(assay=assay, fastq_sets=fastq_sets, metadata=metadata)


class DesignIP(BaseModel):
    """
    Represents an IP (e.g., ChIP/CAT) experiment design:
    paired IP/control FastqSetIP objects, plus per-experiment metadata.

    Attributes:
        assay:       Assay type (e.g., ChIP, CAT).
        experiments: List of IPExperiment, each with .ip and optional .control FastqSetIP.
        metadata:    List of Metadata matching the experiments list.
    """
    assay: Assay
    experiments: list[IPExperiment]
    metadata: list[Metadata]

    @field_validator("assay")
    @classmethod
    def validate_assay(cls, assay: Assay) -> Assay:
        if assay not in Assay.ip_assays():
            raise ValueError(f"Assay '{assay.value}' should use `Design` instead")
        return assay

    @property
    def sample_names(self) -> list[str]:
        """
        All unique IP and control set names, sorted.
        """
        ip_names = [exp.ip_set_fullname for exp in self.experiments]
        ctrl_names = [exp.control_fullname for exp in self.experiments if exp.has_control]
        return sorted(set(ip_names + ctrl_names))

    @property
    def ips_performed(self) -> list[str]:
        """
        Unique IP antibodies/proteins used across experiments.
        """
        return list({exp.ip_performed for exp in self.experiments})

    @property
    def controls_performed(self) -> list[str]:
        """
        Unique control antibodies used across experiments.
        """
        return list({exp.control_performed for exp in self.experiments if exp.has_control})

    def query(
        self,
        sample_name: str,
        full: bool = False,
    ) -> FastqSetIP | dict[str, FastqSetIP]:
        """
        Retrieve IP (and optionally control) FastqSetIP by sample name.

        Args:
            sample_name: IP or control set fullname.
            full: If True, return dict {"ip": ..., "control": ...}; else return the matched set only.

        Raises:
            KeyError: If no match.
        """
        for exp in self.experiments:
            if exp.ip_set_fullname == sample_name:
                out = {"ip": exp.ip, "control": exp.control}
                return out if full else exp.ip
            if exp.has_control and exp.control_fullname == sample_name:
                out = {"ip": exp.ip, "control": exp.control}
                return out if full else exp.control
        raise KeyError(f"Sample '{sample_name}' not found in DesignIP")

    @classmethod
    def from_fastq_files(
        cls,
        assay: Assay,
        files: Iterable[str | pathlib.Path],
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **exp_kwargs: Any,
    ) -> DesignIP:
        """
        Build DesignIP from a mixture of IP/control FASTQ paths.

        Groups by FastqFileIP.sample_base_without_ip and read number.

        Args:
            assay: The assay type (must be an IP assay).
            files: Iterable of file paths (strings or pathlib.Path).
            metadata:
                - Callable(sample_name) → Metadata to customize per-sample metadata.
                - Single Metadata instance applied to all.
                - None → defaults to Metadata(scale_group="all").
            exp_kwargs: Extra fields forwarded to IPExperiment constructor.
        """
        # Convert and sort
        ips: list[FastqFileIP] = [FastqFileIP(path=pathlib.Path(f)) for f in files]

        # Bucket by sample_base_without_ip
        buckets: dict[str, dict[str, list[FastqFileIP]]] = defaultdict(lambda: {"ip": [], "control": []})
        for f in ips:
            key = f.sample_base_without_ip
            side = "control" if f.is_control else "ip"
            buckets[key][side].append(f)

        experiments: list[IPExperiment] = []
        _metadata: list[Metadata] = []
        for name, sides in buckets.items():
            # Sort by read_number
            ip_list = sorted(sides["ip"], key=lambda x: x.read_number)
            ctrl_list = sorted(sides["control"], key=lambda x: x.read_number)

            if not ip_list and not ctrl_list:
                raise ValueError(f"No valid FASTQ files found for sample '{name}'")
            elif not ip_list and ctrl_list:
                ip_list = ctrl_list
                ctrl_list = []
            elif len(ip_list) > 2 or len(ctrl_list) > 2:
                raise ValueError(
                    f"Unexpected number of FASTQ files for '{name}': "
                    f"IP={len(ip_list)}, Control={len(ctrl_list)}"
                )
            
            ip_set = FastqSetIP(name=name, r1=ip_list[0], r2=ip_list[1] if len(ip_list) > 1 else None)
            ctrl_set = (
                FastqSetIP(name=name, r1=ctrl_list[0], r2=ctrl_list[1] if len(ctrl_list) > 1 else None)
                if ctrl_list
                else None
            )

            experiments.append(IPExperiment(ip=ip_set, control=ctrl_set, **exp_kwargs))
            
            # Build Metadata (following the same pattern as Design)
            if callable(metadata):
                meta = metadata(name)
            elif isinstance(metadata, Metadata):
                meta = metadata
            else:
                meta = Metadata(scale_group="all")
            _metadata.append(meta)

        return cls(assay=assay, experiments=experiments, metadata=_metadata)

    @classmethod
    def from_directory(
        cls,
        assay: Assay,
        directory: str | pathlib.Path,
        glob_patterns: Iterable[str] = ("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"),
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **kwargs: Any,
    ) -> DesignIP:  
        """
       Scan a directory for IP/control FASTQ files and build a DesignIP.

        Args:
            directory: Root path to search.
            glob_patterns: Filename patterns to include.
            metadata: Callable(sample_name) → Metadata or single Metadata instance.
            **kwargs: Extra fields converted directly to a shared Metadata.
        """
        dir_path = pathlib.Path(directory)
        files = list(chain.from_iterable(dir_path.glob(p) for p in glob_patterns))
        if not files:
            raise FileNotFoundError(f"No FASTQ files found in {dir_path}")

        if not callable(metadata) and not isinstance(metadata, Metadata):
            # If metadata is not callable, assume it's a single Metadata instance and add the kwargs to it
            metadata = Metadata(**{"scale_group": "all", **kwargs})
        
        return cls.from_fastq_files(assay=assay, files=files, metadata=metadata)

    def to_dataframe(self) -> pd.DataFrame:
        """
        Export IP design to pandas DataFrame, validated by DataFrameDesignIP.

        Columns: sample_name, ip, control, ip_r1, ip_r2, control_r1, control_r2, plus metadata fields.
        """
        import pandas as pd

        rows: list[dict[str, Any]] = []
        for exp, md in zip(self.experiments, self.metadata):
            base = exp.ip.name
            row: dict[str, Any] = {
                "sample_name": base,
                # IP reads
                "r1": exp.ip.r1.path,
                "r2": exp.ip.r2.path if exp.ip.r2 else None,
                # Control reads
                "r1_control": exp.control.r1.path if exp.control else None,
                "r2_control": exp.control.r2.path if exp.control and exp.control.r2 else None,
                # IP/control labels
                "ip": exp.ip_performed,
                "control": exp.control_fullname if exp.has_control else None,
            }
            row.update(md.model_dump(exclude_none=True))
            rows.append(row)

        df = pd.DataFrame(rows).sort_values("sample_name")
        return DataFrame[DesignDataFrame](df)
