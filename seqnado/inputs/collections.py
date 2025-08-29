from __future__ import annotations
import pathlib
from collections import defaultdict
from itertools import chain
from typing import Any, Callable, Iterable
import pandas as pd
from pydantic import BaseModel, field_validator
from pandera.typing import DataFrame
from .fastq import FastqFile, FastqSet, FastqSetIP, FastqFileIP
from .core import Metadata, Assay
from .validation import DesignDataFrame
from .experiment import ExperimentIP


class BaseFastqCollection(BaseModel):
    """
    Base class for all design types providing common functionality.
    """
    assay: Assay
    metadata: list[Metadata]

    @classmethod
    def _build_metadata(
        cls,
        sample_name: str,
        metadata: Callable[[str], Metadata] | Metadata | None,
        assay: Assay,
    ) -> Metadata:
        """Build metadata for a sample ensuring assay is always set.

        Rules:
        - If callable: call with sample_name to get Metadata.
        - If Metadata instance: use directly.
        - If None: create default with norm_group="all".
        In all cases force metadata.assay to the provided assay if it's None or different.
        """
        if callable(metadata):
            md = metadata(sample_name)
        elif isinstance(metadata, Metadata):
            md = metadata
        else:
            md = Metadata(norm_group="all")
        # Always stamp assay (requirement: always include assay in metadata)
        if md.assay != assay:
            md.assay = assay
        return md

    @classmethod
    def _discover_files(
        cls,
        directory: str | pathlib.Path,
        glob_patterns: Iterable[str] = ("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"),
    ) -> list[pathlib.Path]:
        """Discover FASTQ files in a directory using glob patterns."""
        dir_path = pathlib.Path(directory)
        files = list(chain.from_iterable(dir_path.glob(p) for p in glob_patterns))
        if not files:
            raise FileNotFoundError(f"No FASTQ files found in {dir_path}")
        return files

    @classmethod
    def _prepare_metadata_for_directory(
        cls,
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **kwargs: Any,
    ) -> Callable[[str], Metadata] | Metadata:
        """Prepare metadata parameter for directory-based construction."""
        if not callable(metadata) and not isinstance(metadata, Metadata):
            metadata = Metadata(**{"norm_group": "all", **kwargs})
        return metadata
    
    @classmethod
    def from_csv(cls, file_path: str | pathlib.Path) -> FastqCollection:
        """Build a SampleCollection from a CSV file."""
        df = pd.read_csv(file_path)
        return cls.from_dataframe(df)

    def to_dataframe(self) -> pd.DataFrame:
        """Export the design to a pandas DataFrame."""
        raise NotImplementedError("Subclasses must implement to_dataframe")
    
    @property
    def sample_names(self) -> list[str]:
        raise NotImplementedError("Subclasses must implement sample_names")
    
    @property
    def fastq_pairs(self) -> dict[str, list[pathlib.Path]]:
        """Returns a dictionary mapping sample names to their FASTQ file paths."""
        raise NotImplementedError("Subclasses must implement fastq_pairs")

    @property
    def fastq_paths(self) -> list[pathlib.Path]:
        """Flattened list of FASTQ paths from all samples (subclasses must implement get_fastq_files)."""
        return list(chain.from_iterable(self.fastq_pairs.values()))


    def symlink_fastq_files(self, output_dir: str | pathlib.Path) -> None:
        """Symlink FASTQ files to a specified output directory.
        Args:
            output_dir: Directory to create symlinks in.
        """
        output_dir = pathlib.Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        for sample_name, paths in self.fastq_pairs.items():
            for read_number, path in enumerate(paths, start=1):
                symlink_path = output_dir / f"{sample_name}_{read_number}.fastq.gz"
                if not symlink_path.exists():
                    symlink_path.symlink_to(path.resolve().absolute())


class FastqCollection(BaseFastqCollection):
    """
    Represents a collection of sequencing samples (FASTQ files) grouped into named sets,
    with optional per-sample metadata.

    Attributes:
        fastq_sets: List of FastqSet objects (paired or single-end samples).
        metadata:  List of Metadata objects corresponding one-to-one with fastq_sets.
    """
    fastq_sets: list[FastqSet]

    @property
    def primary_file_type(self) -> str:
        return "fastq"

    @field_validator("assay")
    @classmethod
    def validate_non_ip_assay(cls, v: Assay) -> Assay:
        """Ensure the assay doesn't require IP (immunoprecipitation)."""
        if v in Assay.ip_assays():
            raise ValueError(f"Assay '{v.value}' requires IP and should use IPSampleCollection instead")
        return v

    @property
    def sample_ids(self) -> list[str]:
        """
        Returns all sample IDs in the design.
        """
        return [fs.sample_id for fs in self.fastq_sets]
    
    @property
    def sample_names(self) -> list[str]:
        """
        Returns all sample names in the design.
        """
        return self.sample_ids

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

    def get_file_paths(self, kind: str | None = None) -> list[pathlib.Path]:
        kind = kind or "fastq"
        match kind:
            case "fastq":
                return self.fastq_paths
            case "fastq_r1":
                return [fs.r1.path for fs in self.fastq_sets]
            case "fastq_r2":
                return [fs.r2.path for fs in self.fastq_sets if fs.r2]
            case _:
                raise ValueError(f"Unsupported file kind '{kind}' for SampleCollection")
    
    @property
    def fastq_pairs(self) -> dict[str, list[pathlib.Path]]:
        """
        Returns a dictionary mapping sample names to their FASTQ file paths.
        """
        return {
            fs.name: [fs.r1.path, fs.r2.path] if fs.r2 else [fs.r1.path]
            for fs in self.fastq_sets
        }

    def query(self, sample_name: str) -> FastqSet:
        """
        Retrieve the FastqSet by its sample name.

        Raises:
            ValueError: If sample_name not found.
        """
        try:
            return next(fs for fs in self.fastq_sets if fs.sample_id == sample_name)
        except StopIteration:
            raise ValueError(f"Sample '{sample_name}' not found in SampleCollection")

    @classmethod
    def from_fastq_files(
        cls,
        assay: Assay,
        files: Iterable[str | pathlib.Path],
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **fastqset_kwargs: Any,
    ) -> FastqCollection:
        """
        Build a SampleCollection by scanning a list of FASTQ paths:

        1. Convert raw paths to FastqFile.
        2. Group by `sample_base` and sort by read_number.
        3. Create FastqSet (single- or paired-end) for each sample.
        4. Generate Metadata via `metadata(sample_name)`, or default.

        Args:
            files: Iterable of file paths (strings or pathlib.Path).
            metadata:
                - Callable(sample_name) → Metadata to customize per-sample metadata.
                - Single Metadata instance applied to all.
                - None → defaults to Metadata(norm_group="all").
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
                fs = FastqSet(sample_id=sample, r1=fqs[0], **fastqset_kwargs)
            elif len(fqs) == 2:
                fs = FastqSet(sample_id=sample, r1=fqs[0], r2=fqs[1], **fastqset_kwargs)
            else:
                raise ValueError(
                    f"Unexpected number of FASTQ files for '{sample}': {len(fqs)}"
                )
            _fastq_sets.append(fs)

            # Build Metadata using base class method
            _metadata.append(cls._build_metadata(sample, metadata, assay))

        return cls(assay=assay, fastq_sets=_fastq_sets, metadata=_metadata)

    @classmethod
    def from_directory(
        cls,
        assay: Assay,
        directory: str | pathlib.Path,
        glob_patterns: Iterable[str] = ("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"),
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **kwargs: Any,
    ) -> FastqCollection:
        """
        Recursively scan a directory for FASTQ files and build a SampleCollection.

        Args:
            directory: Root path to search.
            glob_patterns: Filename patterns to include.
            metadata: Callable(sample_name) → Metadata or single Metadata instance.
            **kwargs: Extra fields converted directly to a shared Metadata.
        """
        files = cls._discover_files(directory, glob_patterns)
        metadata = cls._prepare_metadata_for_directory(metadata, **kwargs)
        return cls.from_fastq_files(assay=assay, files=files, metadata=metadata)

    def to_dataframe(self, validate: bool = True) -> pd.DataFrame:
        """
        Export the design to a pandas DataFrame, validated by DataFrameDesign.

        Columns: sample_name, r1, r2, plus all metadata fields.
        """
        import pandas as pd

        rows: list[dict[str, Any]] = []
        for fs, md in zip(self.fastq_sets, self.metadata):
            row: dict[str, Any] = {
                "sample_id": fs.sample_id,
                "r1": fs.r1.path,
                "r2": fs.r2.path if fs.r2 else None,
                "uid": f"{fs.sample_id}"
            }
            metadata_dict = md.model_dump(exclude_none=True)
            # Convert Assay enum to string value for schema validation
            if 'assay' in metadata_dict and hasattr(metadata_dict['assay'], 'value'):
                metadata_dict['assay'] = metadata_dict['assay'].value
            row.update(metadata_dict)
            rows.append(row)

        df = pd.DataFrame(rows).sort_values("sample_id").set_index("uid")
        if validate:
            return DataFrame[DesignDataFrame](df)
        else:
            return df

    @classmethod
    def from_dataframe(cls, assay: Assay, df: pd.DataFrame, **fastqset_kwargs: Any) -> FastqCollection:
        """
        Build a SampleCollection from a DataFrame, validated by DataFrameDesign.

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


class FastqCollectionForIP(BaseFastqCollection):
    """
    Represents an IP (e.g., ChIP/CAT) experiment design, consisting of paired IP/control FastqSetIP objects and per-experiment metadata.
        assay (Assay): Assay type (e.g., ChIP, CAT).
        experiments (list[ExperimentIP]): List of IPExperiment objects, each containing .ip and optional .control FastqSetIP.
        metadata (list[Metadata]): List of Metadata objects matching the experiments list.
        sample_names (list[str]): All unique IP and control set names, sorted.
        ips_performed (list[str]): Unique IP antibodies/proteins used across experiments.
        controls_performed (list[str]): Unique control antibodies used across experiments.
        query(sample_name, full=False):
            Scan a directory for IP/control FASTQ files and build an IPSampleCollection.
        to_dataframe():
    Example:
        >>> collection = IPSampleCollection.from_directory(
        ...     assay=Assay.CHIP,
        ...     directory="/data/fastq",
        ...     glob_patterns=["*.fastq.gz"]
        ... )
        >>> df = collection.to_dataframe()
    """
    experiments: list[ExperimentIP]

    @property
    def primary_file_type(self) -> str:
        return "fastq"

    @field_validator("assay")
    @classmethod
    def validate_assay(cls, assay: Assay) -> Assay:
        if assay not in Assay.ip_assays():
            raise ValueError(f"Assay '{assay.value}' should use `SampleCollection` instead")
        return assay

    @property
    def sample_ids(self) -> list[str]:
        """
        All unique IP and control set IDs, sorted.
        """
        ip_names = [exp.ip_set_fullname for exp in self.experiments]
        ctrl_names = [exp.control_fullname for exp in self.experiments if exp.has_control]
        return sorted(set(ip_names + ctrl_names))
    
    @property
    def sample_names(self) -> list[str]:
        """
        All unique IP and control set names, sorted.
        """
        return self.sample_ids

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
    
    @property
    def fastq_pairs(self) -> dict[str, list[pathlib.Path]]:
        """
        Returns a dictionary mapping sample names to their FASTQ file paths.
        """
        files = dict()
        for exp in self.experiments:
            files[exp.ip_set_fullname] = [exp.ip.r1.path, exp.ip.r2.path] if exp.ip.r2 else [exp.ip.r1.path]
            if exp.has_control:
                files[exp.control_fullname] = [exp.control.r1.path, exp.control.r2.path] if exp.control.r2 else [exp.control.r1.path]
        return files

    def get_file_paths(self, kind: str | None = None) -> list[pathlib.Path]:
        kind = kind or "fastq"
        all_fastqs = list({p for paths in self.fastq_pairs.values() for p in paths})
        match kind:
            case "fastq":
                return sorted(all_fastqs)
            case "fastq_r1":
                return sorted({exp.ip.r1.path for exp in self.experiments} | {exp.control.r1.path for exp in self.experiments if exp.control})
            case "fastq_r2":
                return sorted(
                    {p for exp in self.experiments for p in [
                        exp.ip.r2.path if exp.ip.r2 else None,
                        exp.control.r2.path if exp.control and exp.control.r2 else None,
                    ] if p is not None}
                )
            case _:
                raise ValueError(f"Unsupported file kind '{kind}' for SampleCollectionForIP")
    

    def get_control_performed(self, uid: str) -> str:
        """
        Get the control name for the IP sample.
        """
        return self.to_dataframe().loc[uid, "control"]

    def is_paired_end(self, uid: str) -> bool:
        """
        Check if the given sample ID is paired-end.
        """
        return self.to_dataframe().loc[uid, "r2"] is not None

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
        raise KeyError(f"Sample '{sample_name}' not found in IPSampleCollection")

    @classmethod
    def from_fastq_files(
        cls,
        assay: Assay,
        files: Iterable[str | pathlib.Path],
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **exp_kwargs: Any,
    ) -> FastqCollectionForIP:
        """
        Build IPSampleCollection from a mixture of IP/control FASTQ paths.

        Groups by FastqFileIP.sample_base_without_ip and read number.

        Args:
            assay: The assay type (must be an IP assay).
            files: Iterable of file paths (strings or pathlib.Path).
            metadata:
                - Callable(sample_name) → Metadata to customize per-sample metadata.
                - Single Metadata instance applied to all.
                - None → defaults to Metadata(norm_group="all").
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

        experiments: list[ExperimentIP] = []
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
            
            ip_set = FastqSetIP(sample_id=name, r1=ip_list[0], r2=ip_list[1] if len(ip_list) > 1 else None)
            ctrl_set = (
                FastqSetIP(sample_id=name, r1=ctrl_list[0], r2=ctrl_list[1] if len(ctrl_list) > 1 else None)
                if ctrl_list
                else None
            )

            experiments.append(ExperimentIP(ip=ip_set, control=ctrl_set, **exp_kwargs))
            
            # Build Metadata using base class method
            _metadata.append(cls._build_metadata(name, metadata, assay))

        return cls(assay=assay, experiments=experiments, metadata=_metadata)

    @classmethod
    def from_directory(
        cls,
        assay: Assay,
        directory: str | pathlib.Path,
        glob_patterns: Iterable[str] = ("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"),
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **kwargs: Any,
    ) -> FastqCollectionForIP:  
        """
        Scan a directory for IP/control FASTQ files and build a IPSampleCollection.

        Args:
            directory: Root path to search.
            glob_patterns: Filename patterns to include.
            metadata: Callable(sample_name) → Metadata or single Metadata instance.
            **kwargs: Extra fields converted directly to a shared Metadata.
        """
        files = cls._discover_files(directory, glob_patterns)
        metadata = cls._prepare_metadata_for_directory(metadata, **kwargs)
        return cls.from_fastq_files(assay=assay, files=files, metadata=metadata)

    def to_dataframe(self) -> pd.DataFrame:
        """
        Export IP design to pandas DataFrame, validated by DataFrameDesignIP.

        Columns: sample_name, ip, control, ip_r1, ip_r2, control_r1, control_r2, plus metadata fields.
        """
        import pandas as pd

        rows: list[dict[str, Any]] = []
        for exp, md in zip(self.experiments, self.metadata):
            base = exp.ip.sample_id
            row: dict[str, Any] = {
                "sample_id": base,
                # IP reads
                "r1": exp.ip.r1.path,
                "r2": exp.ip.r2.path if exp.ip.r2 else None,
                # Control reads
                "r1_control": exp.control.r1.path if exp.control else None,
                "r2_control": exp.control.r2.path if exp.control and exp.control.r2 else None,
                # IP/control labels
                "ip": exp.ip_performed,
                "control": exp.control_fullname if exp.has_control else None,
                'uid': f"{base}_{exp.ip_performed}"
            }
            row.update(md.model_dump(exclude_none=True, mode="json"))
            rows.append(row)

        df = pd.DataFrame(rows).sort_values("sample_id").set_index('uid')
        return DataFrame[DesignDataFrame](df)
