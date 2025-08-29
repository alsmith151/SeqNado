from __future__ import annotations
import pathlib
from collections import defaultdict
from itertools import chain
from typing import Any, Callable, Iterable, Union
import pandas as pd
from pydantic import BaseModel, field_validator, Field
from pandera.typing import DataFrame
from .fastq import FastqFile, FastqSet, FastqSetIP, FastqFileIP
from .core import Metadata, Assay
from .validation import DesignDataFrame
from .experiment import ExperimentIP


class BaseSampleCollection(BaseModel):
    """
    Base class for all design types providing common functionality.
    """
    assay: Assay
    metadata: list[Metadata]

    @classmethod
    def _build_metadata(
        cls,
        sample_name: str,
        metadata: Callable[[str], Metadata] | Metadata | None = None,
    ) -> Metadata:
        """Build metadata for a sample given the metadata parameter."""
        if callable(metadata):
            return metadata(sample_name)
        elif isinstance(metadata, Metadata):
            return metadata
        else:
            return Metadata(norm_group="all")

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
    def from_csv(cls, file_path: str | pathlib.Path) -> SampleCollection:
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


class SampleCollection(BaseSampleCollection):
    """
    Represents a collection of sequencing samples (FASTQ files) grouped into named sets,
    with optional per-sample metadata.

    Attributes:
        fastq_sets: List of FastqSet objects (paired or single-end samples).
        metadata:  List of Metadata objects corresponding one-to-one with fastq_sets.
    """
    fastq_sets: list[FastqSet]

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
    def fastq_paths(self) -> list[pathlib.Path]:
        """
        Flattens all R1/R2 file paths into a single list.
        """
        return [
            path
            for fs in self.fastq_sets
            for path in (fs.r1.path, *(fs.r2.path if fs.r2 else []))
        ]
    
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
            return next(fs for fs in self.fastq_sets if fs.name == sample_name)
        except StopIteration:
            raise ValueError(f"Sample '{sample_name}' not found in SampleCollection")

    @classmethod
    def from_fastq_files(
        cls,
        assay: Assay,
        files: Iterable[str | pathlib.Path],
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **fastqset_kwargs: Any,
    ) -> SampleCollection:
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
            _metadata.append(cls._build_metadata(sample, metadata))

        return cls(assay=assay, fastq_sets=_fastq_sets, metadata=_metadata)

    @classmethod
    def from_directory(
        cls,
        assay: Assay,
        directory: str | pathlib.Path,
        glob_patterns: Iterable[str] = ("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"),
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **kwargs: Any,
    ) -> SampleCollection:
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
            }
            metadata_dict = md.model_dump(exclude_none=True)
            # Convert Assay enum to string value for schema validation
            if 'assay' in metadata_dict and hasattr(metadata_dict['assay'], 'value'):
                metadata_dict['assay'] = metadata_dict['assay'].value
            row.update(metadata_dict)
            rows.append(row)

        df = pd.DataFrame(rows).sort_values("sample_id")
        if validate:
            return DataFrame[DesignDataFrame](df)
        else:
            return df

    @classmethod
    def from_dataframe(cls, assay: Assay, df: pd.DataFrame, **fastqset_kwargs: Any) -> SampleCollection:
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


class SampleCollectionForIP(BaseSampleCollection):
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
    ) -> SampleCollectionForIP:
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
            
            ip_set = FastqSetIP(name=name, r1=ip_list[0], r2=ip_list[1] if len(ip_list) > 1 else None)
            ctrl_set = (
                FastqSetIP(name=name, r1=ctrl_list[0], r2=ctrl_list[1] if len(ctrl_list) > 1 else None)
                if ctrl_list
                else None
            )

            experiments.append(ExperimentIP(ip=ip_set, control=ctrl_set, **exp_kwargs))
            
            # Build Metadata using base class method
            _metadata.append(cls._build_metadata(name, metadata))

        return cls(assay=assay, experiments=experiments, metadata=_metadata)

    @classmethod
    def from_directory(
        cls,
        assay: Assay,
        directory: str | pathlib.Path,
        glob_patterns: Iterable[str] = ("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"),
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **kwargs: Any,
    ) -> SampleCollectionForIP:  
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


class MultiAssayCollection(BaseSampleCollection):
    """
    Represents a design containing multiple assay types.
    Can handle both regular and IP-based assays in a single design.
    
    Attributes:
        designs: Dictionary mapping assay types to their respective SampleCollection/IPSampleCollection objects.
        metadata: Global metadata list (can be empty if each design has its own metadata).
    """
    designs: dict[Assay, SampleCollection | SampleCollectionForIP]
    
    def __init__(self, designs: dict[Assay, SampleCollection | SampleCollectionForIP], **kwargs):
        # For multi-assay, we set a default assay but it's not really meaningful
        # The real assay information is in the individual designs
        super().__init__(assay=list(designs.keys())[0] if designs else Assay.RNA, 
                        metadata=[], **kwargs)
        self.designs = designs
    
    @property
    def assays(self) -> list[Assay]:
        """Returns all assay types in this multi-assay design."""
        return list(self.designs.keys())
    
    @property
    def sample_names(self) -> list[str]:
        """Returns all sample names across all assays."""
        all_names = []
        for design in self.designs.values():
            all_names.extend(design.sample_names)
        return sorted(set(all_names))
    
    def get_design(self, assay: Assay) -> SampleCollection | SampleCollectionForIP:
        """Get the design for a specific assay type."""
        if assay not in self.designs:
            raise KeyError(f"Assay '{assay.value}' not found in MultiAssayDesign")
        return self.designs[assay]
    
    def query(self, sample_name: str, assay: Assay | None = None) -> FastqSet | FastqSetIP | dict[Assay, FastqSet | FastqSetIP]:
        """
        Query for a sample across assays.
        
        Args:
            sample_name: Name of the sample to find
            assay: Specific assay to search in (if None, searches all)
            
        Returns:
            If assay is specified: FastqSet/FastqSetIP for that sample in that assay
            If assay is None: Dict mapping assays to FastqSet/FastqSetIP objects
        """
        if assay is not None:
            return self.get_design(assay).query(sample_name)
        
        results = {}
        for assay_type, design in self.designs.items():
            try:
                results[assay_type] = design.query(sample_name)
            except (ValueError, KeyError):
                continue  # Sample not found in this assay
        
        if not results:
            raise ValueError(f"Sample '{sample_name}' not found in any assay")
        
        return results
    
    @classmethod
    def from_designs(cls, designs: dict[Assay, SampleCollection | SampleCollectionForIP]) -> MultiAssayCollection:
        """Create a MultiAssayDesign from a dictionary of existing designs."""
        return cls(designs=designs)
    
    @classmethod
    def from_directories(
        cls,
        assay_directories: dict[Assay, str | pathlib.Path],
        glob_patterns: Iterable[str] = ("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"),
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **kwargs: Any,
    ) -> MultiAssayCollection:
        """
        Create a MultiAssayDesign by scanning multiple directories for different assays.
        
        Args:
            assay_directories: Dictionary mapping assay types to their respective directories
            glob_patterns: Filename patterns to include
            metadata: Metadata function or instance to apply to all assays
            **kwargs: Additional arguments passed to individual design constructors
        """
        designs = {}
        
        for assay, directory in assay_directories.items():
            if assay in Assay.ip_assays():
                design = SampleCollectionForIP.from_directory(
                    assay=assay, 
                    directory=directory, 
                    glob_patterns=glob_patterns,
                    metadata=metadata,
                    **kwargs
                )
            else:
                design = SampleCollection.from_directory(
                    assay=assay,
                    directory=directory,
                    glob_patterns=glob_patterns, 
                    metadata=metadata,
                    **kwargs
                )
            designs[assay] = design
        
        return cls(designs=designs)
    
    def to_dataframe(self) -> pd.DataFrame:
        """
        Export the multi-assay design to a pandas DataFrame.
        Combines all assays into a single DataFrame with an 'assay' column.
        """
        import pandas as pd
        
        all_dfs = []
        for assay, design in self.designs.items():
            df = design.to_dataframe()
            df['assay'] = assay.value
            all_dfs.append(df)
        
        if not all_dfs:
            return pd.DataFrame()
        
        combined_df = pd.concat(all_dfs, ignore_index=True)
        return combined_df.sort_values(['assay', 'sample_name'])
    
    def split_by_assay(self) -> dict[Assay, SampleCollection | SampleCollectionForIP]:
        """Return the individual designs by assay type."""
        return self.designs.copy()



SampleCollectionType = Union[SampleCollection, SampleCollectionForIP, MultiAssayCollection]

def select_sample_collection(
    assay: Assay,
    **kwargs: Any,
) -> SampleCollectionType:
    """
    Select the appropriate SampleCollection type based on the assay.
    Args:
        assay: The assay type (e.g., Assay.RNA, Assay.CHIP).
        **kwargs: Additional parameters passed to the collection constructor.
    
    Returns:
        An instance of SampleCollection, SampleCollectionForIP, or MultiAssayCollection.
    """
    if assay in Assay.ip_assays():
        return SampleCollectionForIP(assay=assay, **kwargs)
    else:
        return SampleCollection(assay=assay, **kwargs)




class SampleGroup(BaseModel):
    """A single group of samples with an optional reference sample."""
    name: str
    samples: list[str]
    reference_sample: str | None = None

    def __len__(self): return len(self.samples)
    def __contains__(self, sample): return sample in self.samples
    def __str__(self): return f"{self.name}: {len(self.samples)} samples (ref={self.reference_sample})"


class SampleGroups(BaseModel):
    """
    A collection of SampleGroup instances that represent one grouping scheme
    (e.g., for normalization or scaling).
    """
    groups: list[SampleGroup]

    @classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        subset_column: str = "norm_group",
        *,
        reference_sample: str | None = None,
        include_controls: bool = False,
    ) -> "SampleGroups":
        """
        Build multiple SampleGroups from a DataFrame based on a grouping column.
        """
        if subset_column not in df.columns:
            raise ValueError(f"Column '{subset_column}' not found in design.")

        groups = []
        for group_value, group_df in df.groupby(subset_column):
            sample_names = group_df.index.tolist()
            ref_sample = reference_sample if reference_sample in sample_names else sample_names[0]
            groups.append(SampleGroup(name=str(group_value), samples=sample_names, reference_sample=ref_sample))

        return cls(groups=groups)

    def sample_to_group(self) -> dict[str, str]:
        return {sample: g.name for g in self.groups for sample in g.samples}

    def group_to_samples(self) -> dict[str, list[str]]:
        return {g.name: g.samples for g in self.groups}

    def get_group(self, name: str) -> SampleGroup:
        for group in self.groups:
            if group.name == name:
                return group
        raise KeyError(f"Group '{name}' not found.")

    def get_samples(self, group_name: str) -> list[str]:
        return self.get_group(group_name).samples
    
    @property
    def empty(self) -> bool:
        """Check if there are no groups defined."""
        return len(self.groups) == 0

    def __str__(self):
        return f"SampleGroups({len(self.groups)} groups: {', '.join(g.name for g in self.groups)})"
    
    def __len__(self):
        return len(self.groups)


class SampleGroupings(BaseModel):
    """
    A container for multiple named SampleGroups sets,
    e.g., for 'normalization', 'scaling', 'visualization', etc.
    """
    groupings: dict[str, SampleGroups] = Field(default_factory=dict)

    def add_grouping(self, name: str, groups: SampleGroups):
        self.groupings[name] = groups

    def get_grouping(self, name: str) -> SampleGroups:
        if name not in self.groupings:
            raise KeyError(f"Grouping '{name}' not found.")
        return self.groupings[name]

    def __contains__(self, name: str) -> bool:
        return name in self.groupings

    def __str__(self):
        return f"SampleGroupings({', '.join(self.groupings)})"

    def all_samples(self) -> list[str]:
        """Return all unique samples across all groupings."""
        return list({sample for g in self.groupings.values() for group in g.groups for sample in group.samples})
    
    @property
    def empty(self) -> bool:
        """Check if there are no groupings defined."""
        return len(self.groupings) == 0

