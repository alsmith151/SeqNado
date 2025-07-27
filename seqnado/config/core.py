from pydantic import BaseModel, field_validator, computed_field
from datetime import datetime
from typing import Union, Literal
from pathlib import Path
from enum import Enum
from seqnado.inputs import Assay
from seqnado import (
    PileupMethod,
    PCRDuplicateTool,
    PCRDuplicateHandling,
    PeakCallingMethod,
    RNAQuantificationMethod,
    SNPCallingMethod,
    MethylationMethod,
    SpikeInMethod,
)


class UserFriendlyError(Exception):
    """Custom exception for user-friendly error messages that should be displayed without traceback."""

    pass


class BowtieIndex(BaseModel):
    """Container for Bowtie2 index files."""

    type: Literal["Bowtie2"] = "Bowtie2"
    prefix: str

    @field_validator("prefix")
    def validate_prefix(cls, v: str) -> str:
        # The prefix should be a valid path with a prefix for the index files.
        # Ensure that its parent is a valid directory.
        # Ensure that index files exist when the prefix is set.
        p = Path(v)
        if not p.parent.exists() and not p.parent.is_dir():
            raise ValueError(
                f"The directory {p.parent} does not exist. Please provide a valid path together with the prefix."
            )

        # Check if index files exist
        for suffix in [
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ]:
            if not (p.with_suffix(suffix).exists() and p.with_suffix(suffix).is_file()):
                raise ValueError(
                    f"Index file {p.with_suffix(suffix)} does not exist. Please provide a valid path together with the prefix."
                )

        return v

    @property
    def files(self) -> list[Path]:
        """Return a list of all Bowtie2 index files."""
        return [
            Path(self.prefix + suffix)
            for suffix in [
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ]
        ]


class STARIndex(BaseModel):
    """Container for STAR index files."""

    type: Literal["STAR"] = "STAR"
    prefix: Path

    @field_validator("prefix")
    def validate_prefix(cls, v: Path) -> Path:
        # The prefix should be a valid path with a prefix for the index files.
        # Ensure that its parent is a valid directory.

        if not v.exists() or not v.is_dir():
            raise ValueError(
                f"The directory {v} does not exist or is not a directory. Please provide a valid path for the STAR index prefix."
            )

    @property
    def files(self) -> list[Path]:
        """Return a list of all STAR index files."""
        return list(self.prefix.glob("*"))


class GenomeConfig(BaseModel):
    """Configuration for genome-related files and indices."""

    name: str
    index: BowtieIndex | STARIndex
    fasta: Path | None = None
    chromosome_sizes: Path | None = None
    gtf: Path | None = None
    genes: Path | None = None
    blacklist: Path | None = None

    organism: str | None = None
    version: str | None = None

    @field_validator("fasta", "chromosome_sizes", "gtf", "genes", "blacklist")
    def validate_paths_exist(cls, v: Path | None) -> Path | None:
        if v is not None and not v.exists():
            # Get the field name from validation info if available
            import inspect

            frame = inspect.currentframe()
            field_name = "File"
            try:
                # Try to get the field name from the call stack
                caller_locals = frame.f_back.f_locals
                if "field_name" in caller_locals:
                    field_name = caller_locals["field_name"].title()
            except Exception as e:
                pass
            finally:
                del frame

            # Handle both relative and absolute paths properly
            original_path = str(v)
            try:
                abs_path = v.resolve()
                parent_dir = abs_path.parent
            except:
                # Fallback for problematic paths
                abs_path = v
                parent_dir = v.parent

            error_msg = "❌ Genome file not found\n\n"
            error_msg += f"📁 Looking for: {v.name}\n"
            error_msg += f"📂 In directory: {parent_dir}\n"
            error_msg += f"� Full path: {original_path}\n\n"

            if parent_dir.exists():
                error_msg += "💡 Directory exists, but file is missing\n\n"

                # Show available files in the directory
                try:
                    all_files = [f.name for f in parent_dir.iterdir() if f.is_file()]
                    if all_files:
                        # Show files with similar names first
                        filename_lower = v.name.lower()
                        similar_files = [
                            f
                            for f in all_files
                            if filename_lower in f.lower()
                            or f.lower() in filename_lower
                        ]
                        other_files = [f for f in all_files if f not in similar_files]

                        if similar_files:
                            error_msg += "🎯 Similar files found:\n"
                            for file in sorted(similar_files)[:5]:
                                error_msg += f"   • {file}\n"
                            error_msg += "\n"

                        if other_files and len(all_files) <= 15:
                            error_msg += "📄 All files in directory:\n"
                            for file in sorted(other_files)[:10]:
                                error_msg += f"   • {file}\n"
                            if len(other_files) > 10:
                                error_msg += (
                                    f"   ... and {len(other_files) - 10} more files\n"
                                )
                        elif other_files:
                            error_msg += f"📄 Directory contains {len(other_files)} other files\n"
                        error_msg += "\n"
                    else:
                        error_msg += "� Directory is empty\n\n"
                except Exception:
                    error_msg += "📄 Could not list directory contents\n\n"
            else:
                error_msg += "💡 Directory does not exist\n\n"

            error_msg += "🔧 How to fix:\n"
            error_msg += "   1. Check the file path for typos\n"
            error_msg += "   2. Ensure the file has been created/downloaded\n"
            error_msg += "   3. Verify file permissions\n"
            error_msg += "   4. Use absolute paths to avoid confusion\n"
            error_msg += f"   5. Current working directory: {Path.cwd()}"

            raise UserFriendlyError(error_msg)
        return v


class ProjectConfig(BaseModel):
    """Configuration for the SeqNado project."""

    name: str
    date: datetime
    directory: Path


class PCRDuplicatesConfig(BaseModel):
    """Configuration for PCR duplicates handling."""

    strategy: PCRDuplicateHandling = PCRDuplicateHandling.NONE
    tool: PCRDuplicateTool | None = None


class QCConfig(BaseModel):
    """Configuration for library quality control."""

    calculate_library_size: bool = True
    calculate_fraction_of_reads_in_peaks: bool = True


class ComputedFieldsMixin:
    """Mixin class providing common computed fields for assay configurations."""

    @computed_field
    @property
    def create_bigwigs(self) -> bool:
        """Whether to make bigwigs (computed from bigwigs config presence)."""
        return getattr(self, "bigwigs", None) is not None

    @computed_field
    @property
    def plot_with_plotnado(self) -> bool:
        """Whether to perform plotting (computed from plotting config presence)."""
        return getattr(self, "plotting", None) is not None

    @computed_field
    @property
    def create_dataset(self) -> bool:
        """Whether to make dataset."""
        return getattr(self, "dataset_for_ml", None) is not None

    @computed_field
    @property
    def create_ucsc_hub(self) -> bool:
        """Whether to make UCSC hub (computed from ucsc_hub config presence)."""
        return getattr(self, "ucsc_hub", None) is not None


class PeakCallingMixin:
    """Mixin for assays that support peak calling."""

    @computed_field
    @property
    def call_peaks(self) -> bool:
        """Whether to call peaks (computed from peak_calling config presence)."""
        return getattr(self, "peak_calling", None) is not None


class SNPCallingMixin:
    """Mixin for assays that support SNP calling."""

    @computed_field
    @property
    def call_snps(self) -> bool:
        """Whether to call SNPs (computed from snp_calling config presence)."""
        return getattr(self, "snp_calling", None) is not None


class MethylationMixin:
    """Mixin for assays that support methylation calling."""

    @computed_field
    @property
    def call_methylation(self) -> bool:
        """Whether to call methylation (computed from methylation config presence)."""
        return getattr(self, "methylation", None) is not None


class PathValidatorMixin:
    """Mixin providing common path validation methods."""

    @staticmethod
    def validate_path_exists(
        v: Path | str | None, field_name: str = "path"
    ) -> Path | str | None:
        """Validate that a path exists if provided."""
        if v is not None:
            path = Path(v) if isinstance(v, str) else v
            if not path.exists():
                raise ValueError(f"{field_name} {v} does not exist.")
        return v

    @staticmethod
    def validate_required_when(
        v: str | None, condition: bool, field_name: str
    ) -> str | None:
        """Validate that a field is provided when a condition is met."""
        if condition and not v:
            raise ValueError(f"{field_name} must be provided when condition is met.")
        return v


class BigwigConfig(BaseModel):
    """Configuration for bigwig generation."""

    pileup_method: list[PileupMethod] | None = None
    binsize: int | None = None


class PlottingConfig(BaseModel):
    """Configuration for plotting features."""

    coordinates: str | None = None
    genes: str | None = None


class PeakCallingConfig(BaseModel):
    """Configuration for peak calling."""

    method: list[PeakCallingMethod] | None = None
    consensus_counts: bool = False


class SpikeInConfig(BaseModel):
    """Configuration for spike-in normalization."""

    method: SpikeInMethod
    exogenous_genome: str | None = None
    endogenous_genome: str | None = None


class UCSCHubConfig(BaseModel):
    """Configuration for UCSC Hub generation."""

    directory: str = "seqnado_output/hub/"
    email: str
    color_by: str = "samplename"


class RNAQuantificationConfig(BaseModel, PathValidatorMixin):
    """Configuration for RNA quantification."""

    method: RNAQuantificationMethod
    salmon_index: str | None = None
    run_deseq2: bool = False

    @field_validator("salmon_index")
    def validate_salmon_index(cls, v: str | None) -> str | None:
        return cls.validate_path_exists(v, "Salmon index path")


class SNPCallingConfig(BaseModel, PathValidatorMixin):
    """Configuration for SNP calling."""

    method: SNPCallingMethod
    annotate_snps: bool = False
    snp_database: str | None = None

    @field_validator("snp_database")
    def validate_snp_database(cls, v: str | None, info) -> str | None:
        annotate_snps = info.data.get("annotate_snps", False)
        return cls.validate_required_when(v, annotate_snps, "snp_database")


class MCCConfig(BaseModel, PathValidatorMixin):
    """Configuration for MCC (Capture-C) analysis."""

    viewpoints: Path
    resolutions: list[int] = [100]

    @field_validator("viewpoints")
    def validate_viewpoints(cls, v: Path) -> Path:
        return cls.validate_path_exists(v, "Viewpoints path")

    @field_validator("resolutions")
    def validate_resolutions(cls, v: list[int]) -> list[int]:
        if not v:
            raise ValueError("Resolutions list must not be empty.")
        if any(resolution <= 0 for resolution in v):
            raise ValueError("All resolutions must be positive integers.")
        return v


class MethylationConfig(BaseModel):
    """Configuration for methylation calling."""

    assay: MethylationMethod


class MLDatasetConfig(BaseModel, PathValidatorMixin):
    """Configuration for ML dataset generation."""

    regions_bed: Path | None = None
    binsize: int | None = None

    @field_validator("regions_bed")
    def validate_regions_bed(cls, v: Path | None) -> Path | None:
        return cls.validate_path_exists(v, "Regions BED file")

    @field_validator("binsize")
    def validate_binsize(cls, v: int | None) -> int | None:
        if v is not None and v <= 0:
            raise ValueError("Binsize must be a positive integer.")
        return v

    def model_post_init(self, __context) -> None:
        """Validate that at least one of regions_bed or binsize is provided."""
        if self.regions_bed is None and self.binsize is None:
            raise ValueError("At least one of regions_bed or binsize must be provided.")


class BaseAssayConfig(BaseModel, ComputedFieldsMixin):
    """Base configuration for all assays."""

    bigwigs: BigwigConfig | None = None
    plotting: PlottingConfig | None = None
    ucsc_hub: UCSCHubConfig | None = None
    dataset_for_ml: MLDatasetConfig | None = None

    # Boolean flags for optional features
    create_geo_submission_files: bool = False


class ATACAssayConfig(BaseAssayConfig, PeakCallingMixin):
    """Configuration specific to ATAC-seq assays."""

    tn5_shift: bool = False
    peak_calling: PeakCallingConfig | None = None
    create_heatmaps: bool = False


class ChIPAssayConfig(BaseAssayConfig, PeakCallingMixin):
    """Configuration specific to ChIP-seq assays."""

    spikein: SpikeInConfig | None = None
    peak_calling: PeakCallingConfig | None = None
    create_heatmaps: bool = False


class CATAssayConfig(BaseAssayConfig, PeakCallingMixin):
    """Configuration specific to CAT-seq assays."""

    spikein: SpikeInConfig | None = None
    tn5_shift: bool = False
    peak_calling: PeakCallingConfig | None = None
    create_heatmaps: bool = False


class RNAAssayConfig(BaseAssayConfig):
    """Configuration specific to RNA-seq assays."""

    rna_quantification: RNAQuantificationConfig | None = None
    create_heatmaps: bool = False


class SNPAssayConfig(BaseAssayConfig, SNPCallingMixin):
    """Configuration specific to SNP calling assays."""

    snp_calling: SNPCallingConfig | None = None


class MCCAssayConfig(BaseAssayConfig):
    """Configuration specific to MCC (Capture-C) assays."""

    mcc: MCCConfig | None = None


class METHAssayConfig(BaseAssayConfig, MethylationMixin):
    """Configuration specific to methylation assays."""

    methylation: MethylationConfig | None = None


class CRISPRAssayConfig(BaseAssayConfig):
    """Configuration specific to CRISPR assays."""

    # CRISPR-specific options can be added here
    pass


# Union type for all assay-specific configurations
AssaySpecificConfig = Union[
    ATACAssayConfig,
    ChIPAssayConfig,
    CATAssayConfig,
    RNAAssayConfig,
    SNPAssayConfig,
    MCCAssayConfig,
    METHAssayConfig,
    CRISPRAssayConfig,
]

# Class constant for assay config mapping to reduce redundancy
ASSAY_CONFIG_MAP = {
    Assay.ATAC: ATACAssayConfig,
    Assay.CHIP: ChIPAssayConfig,
    Assay.CAT: CATAssayConfig,
    Assay.RNA: RNAAssayConfig,
    Assay.SNP: SNPAssayConfig,
    Assay.MCC: MCCAssayConfig,
    Assay.METH: METHAssayConfig,
    Assay.CRISPR: CRISPRAssayConfig,
}


class AssayConfig(Enum):
    ATAC = ATACAssayConfig
    CHIP = ChIPAssayConfig
    CAT = CATAssayConfig
    RNA = RNAAssayConfig
    SNP = SNPAssayConfig
    MCC = MCCAssayConfig
    METH = METHAssayConfig
    CRISPR = CRISPRAssayConfig


class WorkflowConfig(BaseModel):
    """Configuration for the SeqNado workflow."""

    assay: Assay
    project: ProjectConfig
    genome: GenomeConfig
    metadata: Path

    pcr_duplicates: PCRDuplicatesConfig = PCRDuplicatesConfig()
    qc: QCConfig = QCConfig()
    assay_config: AssaySpecificConfig | None = None
    options: str | None = None

    @field_validator("assay_config", mode="before")
    def validate_assay_config_matches_assay(cls, v, info):
        """Ensure the assay_config type matches the specified assay."""
        if v is None:
            return v

        assay = info.data.get("assay")
        if assay is None:
            return v

        expected_config_class = ASSAY_CONFIG_MAP.get(assay)
        if expected_config_class and not isinstance(v, expected_config_class):
            if isinstance(v, dict):
                # Try to create the appropriate config from dict
                return expected_config_class(**v)
            else:
                raise ValueError(
                    f"assay_config must be of type {expected_config_class.__name__} for assay {assay.value}"
                )

        return v

    @classmethod
    def create_assay_config(cls, assay: Assay, **kwargs) -> AssaySpecificConfig:
        """Create the appropriate assay config for the given assay type."""
        config_class = cls._ASSAY_CONFIG_MAP.get(assay)
        if config_class is None:
            raise ValueError(
                f"No configuration class available for assay {assay.value}"
            )

        return config_class(**kwargs)
