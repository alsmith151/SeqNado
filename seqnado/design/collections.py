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

    @property
    def sample_names(self) -> list[str]:
        """Returns all sample names in the design."""
        raise NotImplementedError("Subclasses must implement sample_names")

    def to_dataframe(self) -> pd.DataFrame:
        """Export the design to a pandas DataFrame."""
        raise NotImplementedError("Subclasses must implement to_dataframe")


class SampleCollection(BaseFastqCollection):
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
                fs = FastqSet(name=sample, r1=fqs[0], **fastqset_kwargs)
            elif len(fqs) == 2:
                fs = FastqSet(name=sample, r1=fqs[0], r2=fqs[1], **fastqset_kwargs)
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

    def to_dataframe(self) -> pd.DataFrame:
        """
        Export the design to a pandas DataFrame, validated by DataFrameDesign.

        Columns: sample_name, r1, r2, plus all metadata fields.
        """
        import pandas as pd

        rows: list[dict[str, Any]] = []
        for fs, md in zip(self.fastq_sets, self.metadata):
            row: dict[str, Any] = {
                "sample_name": fs.sample_id,
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


class IPSampleCollection(BaseFastqCollection):
    """
    Represents an IP (e.g., ChIP/CAT) experiment design:
    paired IP/control FastqSetIP objects, plus per-experiment metadata.

    Attributes:
        assay:       Assay type (e.g., ChIP, CAT).
        experiments: List of IPExperiment, each with .ip and optional .control FastqSetIP.
        metadata:    List of Metadata matching the experiments list.
    """
    experiments: list[ExperimentIP]

    @field_validator("assay")
    @classmethod
    def validate_assay(cls, assay: Assay) -> Assay:
        if assay not in Assay.ip_assays():
            raise ValueError(f"Assay '{assay.value}' should use `SampleCollection` instead")
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
        raise KeyError(f"Sample '{sample_name}' not found in IPSampleCollection")

    @classmethod
    def from_fastq_files(
        cls,
        assay: Assay,
        files: Iterable[str | pathlib.Path],
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **exp_kwargs: Any,
    ) -> IPSampleCollection:
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
    ) -> IPSampleCollection:  
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


class MultiAssayDesign(BaseFastqCollection):
    """
    Represents a design containing multiple assay types.
    Can handle both regular and IP-based assays in a single design.
    
    Attributes:
        designs: Dictionary mapping assay types to their respective SampleCollection/IPSampleCollection objects.
        metadata: Global metadata list (can be empty if each design has its own metadata).
    """
    designs: dict[Assay, SampleCollection | IPSampleCollection]
    
    def __init__(self, designs: dict[Assay, SampleCollection | IPSampleCollection], **kwargs):
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
    
    def get_design(self, assay: Assay) -> SampleCollection | IPSampleCollection:
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
    def from_designs(cls, designs: dict[Assay, SampleCollection | IPSampleCollection]) -> MultiAssayDesign:
        """Create a MultiAssayDesign from a dictionary of existing designs."""
        return cls(designs=designs)
    
    @classmethod
    def from_directories(
        cls,
        assay_directories: dict[Assay, str | pathlib.Path],
        glob_patterns: Iterable[str] = ("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"),
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **kwargs: Any,
    ) -> MultiAssayDesign:
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
                design = IPSampleCollection.from_directory(
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
    
    def split_by_assay(self) -> dict[Assay, SampleCollection | IPSampleCollection]:
        """Return the individual designs by assay type."""
        return self.designs.copy()




class SampleGroup(BaseModel):
    """
    Represents a normalization group containing samples and an optional reference sample.
    
    Attributes:
        group: Group identifier (name or number), defaults to "all"
        samples: List of sample names in this group
        reference_sample: Optional reference sample for normalization
    """
    group: str | int = "all"
    samples: list[str]
    reference_sample: str | None = None

    @classmethod
    def from_design(
        cls,
        design: SampleCollection | IPSampleCollection,
        reference_sample: str | None = None,
        subset_column: str = "norm_group",
        subset_value: list[str] | None = None,
        include_controls: bool = False,
    ) -> SampleGroup:
        """
        Create a SampleGroup from a design object.
        
        Args:
            design: The design object (SampleCollection or IPSampleCollection)
            reference_sample: Optional reference sample name
            subset_column: Column to use for subsetting (default: "norm_group")
            subset_value: Values to filter by in the subset column
            include_controls: Whether to include control samples (for DesignIP)
            
        Returns:
            SampleGroup instance
            
        Raises:
            ValueError: If no samples found after filtering
        """
        df = cls._prepare_dataframe(design, include_controls)
        
        if subset_value:
            # Check if subset_column exists in the DataFrame
            if subset_column not in df.columns:
                raise ValueError(f"No samples found for group with {subset_column} in {subset_value}")
            
            df = df.query(f"{subset_column} in {subset_value}")
        
        samples = df.index.tolist()
        if not samples:
            raise ValueError(f"No samples found for group with {subset_column} in {subset_value}")
        
        final_reference = reference_sample or samples[0]
        group_name = subset_value[0] if subset_value else "all"
        
        return cls(
            samples=samples,
            reference_sample=final_reference,
            group=group_name,
        )
    
    @staticmethod
    def _prepare_dataframe(design: SampleCollection | IPSampleCollection, include_controls: bool) -> pd.DataFrame:
        """
        Prepare a DataFrame from design with appropriate sample naming.
        
        Args:
            design: The design object
            include_controls: Whether to include control samples
            
        Returns:
            DataFrame with sample_fullname as index
        """
        base_df = design.to_dataframe()
        
        if isinstance(design, SampleCollection):
            return (
                base_df
                .assign(sample_fullname=lambda df: df.sample_name)
                .set_index("sample_fullname")
            )
        
        # Handle IPSampleCollection
        if not include_controls:
            return (
                base_df
                .assign(sample_fullname=lambda df: df.sample_name + "_" + df.ip)
                .set_index("sample_fullname")
            )
        
        # Include both IP and control samples
        df_ip = (
            base_df
            .assign(sample_fullname=lambda df: df.sample_name + "_" + df.ip)
            .set_index("sample_fullname")
        )
        df_control = (
            base_df
            .query("control.notnull()")
            .assign(sample_fullname=lambda df: df.sample_name + "_" + df.control)
            .set_index("sample_fullname")
        )
        return pd.concat([df_ip, df_control])
    
    def __len__(self) -> int:
        """Return the number of samples in this group."""
        return len(self.samples)
    
    def __contains__(self, sample: str) -> bool:
        """Check if a sample is in this group."""
        return sample in self.samples
    
    def __str__(self) -> str:
        """String representation of the group."""
        return f"SampleGroup(group='{self.group}', samples={len(self.samples)}, reference='{self.reference_sample}')"


class SampleGroups(BaseModel):
    """
    Collection of normalization groups with utilities for sample group management.
    
    Attributes:
        groups: List of SampleGroup objects
    """
    groups: list[SampleGroup]

    @classmethod
    def from_design(
        cls,
        design: SampleCollection | IPSampleCollection,
        reference_sample: str | None = None,
        subset_column: str = "norm_group",
        include_controls: bool = False,
    ) -> SampleGroups:
        """
        Create SampleGroups from a design object.
        
        Args:
            design: The design object (SampleCollection or IPSampleCollection)
            reference_sample: Optional reference sample name
            subset_column: Column to use for grouping (default: "norm_group")
            include_controls: Whether to include control samples (for DesignIP)
            
        Returns:
            SampleGroups instance
        """
        df = design.to_dataframe()

        # Create groups based on subset column if it exists
        if subset_column in df.columns:
            unique_values = df[subset_column].dropna().unique()
            
            groups = [
                SampleGroup.from_design(
                    design=design,
                    reference_sample=reference_sample,
                    subset_column=subset_column,
                    subset_value=[str(value)],
                    include_controls=include_controls,
                )
                for value in unique_values
            ]
        else:
            # Create single group with all samples
            groups = [
                SampleGroup(
                    group="all",
                    samples=design.sample_names,
                    reference_sample=reference_sample or design.sample_names[0],
                )
            ]

        return cls(groups=groups)

    @property
    def sample_to_group_mapping(self) -> dict[str, str | int]:
        """
        Map each sample to its group identifier.
        
        Returns:
            Dictionary mapping sample names to group identifiers
        """
        return {
            sample: group.group 
            for group in self.groups 
            for sample in group.samples
        }

    @property
    def group_to_samples_mapping(self) -> dict[str | int, list[str]]:
        """
        Map each group identifier to its samples.
        
        Returns:
            Dictionary mapping group identifiers to lists of sample names
        """
        return {group.group: group.samples for group in self.groups}

    def get_sample_group(self, sample: str) -> str | int:
        """
        Get the group identifier for a specific sample.
        
        Args:
            sample: Sample name to look up
            
        Returns:
            Group identifier
            
        Raises:
            KeyError: If sample not found in any group
        """
        mapping = self.sample_to_group_mapping
        if sample not in mapping:
            raise KeyError(f"Sample '{sample}' not found in any normalization group")
        return mapping[sample]

    def get_samples_in_group(self, group: str | int) -> list[str]:
        """
        Get all samples in a specific group.
        
        Args:
            group: Group identifier
            
        Returns:
            List of sample names in the group
            
        Raises:
            KeyError: If group not found
        """
        mapping = self.group_to_samples_mapping
        if group not in mapping:
            available_groups = list(mapping.keys())
            raise KeyError(f"Group '{group}' not found. Available groups: {available_groups}")
        return mapping[group]
    
    def get_group_by_name(self, group_name: str | int) -> SampleGroup:
        """
        Get a SampleGroup object by its group identifier.
        
        Args:
            group_name: Group identifier to find
            
        Returns:
            SampleGroup object
            
        Raises:
            KeyError: If group not found
        """
        for group in self.groups:
            if group.group == group_name:
                return group
        raise KeyError(f"Group '{group_name}' not found")
    
    def __len__(self) -> int:
        """Return the number of groups."""
        return len(self.groups)
    
    def __iter__(self):
        """Iterate over groups."""
        return iter(self.groups)
    
    def __contains__(self, group_name: str | int) -> bool:
        """Check if a group exists."""
        return group_name in self.group_to_samples_mapping
    
    def __str__(self) -> str:
        """String representation of the groups."""
        group_info = [f"'{g.group}': {len(g.samples)} samples" for g in self.groups]
        return f"SampleGroups({len(self.groups)} groups: {', '.join(group_info)})"

    # Backward compatibility properties
    @property  
    def sample_groups(self) -> dict[str | int, list[str]]:
        """Backward compatibility alias for group_to_samples_mapping."""
        return self.group_to_samples_mapping

    @property
    def group_samples(self) -> dict[str, str | int]:
        """Backward compatibility alias for sample_to_group_mapping."""
        return self.sample_to_group_mapping

    def get_grouped_samples(self, group: str | int) -> list[str]:
        """Backward compatibility alias for get_samples_in_group."""
        return self.get_samples_in_group(group)
