import re
from pydantic import BaseModel, field_validator, computed_field, Field
from datetime import date as _date
from typing import Union, Literal
from pathlib import Path
from enum import Enum
from seqnado import Assay
from seqnado import (
    PileupMethod,
    PCRDuplicateTool,
    PCRDuplicateHandling,
    PeakCallingMethod,
    QuantificationMethod,
    SNPCallingMethod,
    MethylationMethod,
    SpikeInMethod,
)
from .mixins import (
    CommonComputedFieldsMixin,
    PeakCallingMixin,
    SNPCallingMixin,
    MethylationMixin,
    PathValidatorMixin,
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
    bin_size: int | None = None
    organism: str | None = None
    version: str | None = None


    def model_post_init(self, context):
        if not self.organism:
            self.organism = self.predict_organism()

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

    def predict_organism(self) -> str:
        if "hg" in self.name:
            return "Homo sapiens"
        elif "mm" in self.name:
            return "Mus musculus"
        else:
            return "Unknown"

class ProjectConfig(BaseModel):
    """Configuration for the SeqNado project."""

    name: str
    date: _date = Field(default_factory=_date.today)
    description: str | None = None


class PCRDuplicatesConfig(BaseModel):
    """Configuration for PCR duplicates handling."""

    strategy: PCRDuplicateHandling = PCRDuplicateHandling.NONE
    tool: PCRDuplicateTool | None = None


class QCConfig(BaseModel):
    """Configuration for library quality control."""
    run_fastq_screen: bool = True
    calculate_library_size: bool = True
    calculate_fraction_of_reads_in_peaks: bool = True



class BigwigConfig(BaseModel):
    """Configuration for bigwig generation."""

    pileup_method: list[PileupMethod] | None = None
    binsize: int | None = None


class PlottingConfig(BaseModel):
    """Configuration for plotting features."""

    coordinates: str | None = None
    genes: str | None = None
    file_format: Literal["png", "pdf", "svg"] = "pdf"


class PeakCallingConfig(BaseModel):
    """Configuration for peak calling."""

    method: list[PeakCallingMethod] | None = None
    consensus_counts: bool = False


class SpikeInConfig(BaseModel):
    """Configuration for spike-in normalization."""

    method: SpikeInMethod
    exogenous_genome: str | None = None
    endogenous_genome: str | None = None


class GenomicCoordinate(BaseModel):
    """Configuration for genomic coordinates."""

    chromosome: str
    start: int
    end: int

    @field_validator("start", "end")
    def validate_coordinates(cls, v: int) -> int:
        if v < 0:
            raise ValueError("Genomic coordinates must be non-negative.")
        return v

    # Check that end is greater than start
    @field_validator("end")
    def validate_end_greater_than_start(cls, v: int, info) -> int:
        if v < info.data["start"]:
            raise ValueError("End coordinate must be greater than start coordinate.")
        return v

    @classmethod
    def from_string(cls, coord_str: str) -> "GenomicCoordinate":
        """
        Create a GenomicCoordinate instance from a string representation.
        """
        chromosome, positions = coord_str.split(":")
        start, end = map(int, positions.split("-"))
        return cls(chromosome=chromosome, start=start, end=end)

class UCSCHubConfig(BaseModel):
    """Configuration for UCSC Hub generation."""

    directory: str = "seqnado_output/hub/"
    name: str = "seqnado_hub"
    genome: str = "hg38"
    email: str = "test@example.com"
    two_bit: str | None = None
    organism: str | None = None
    default_position: GenomicCoordinate = Field(default_factory=GenomicCoordinate.from_string("chr1:100-200"), description="Default genomic position")
    color_by: list[str] = Field(default_factory=lambda: ['samplename'], description="List of fields to color the bigwigs")
    overlay_by: list[str] = Field(default_factory=list, description="List of fields to overlay the bigwigs")
    subgroup_by: list[str] = Field(default_factory=lambda: ['method', "norm"], description="List of fields to subgroup the bigwigs")
    supergroup_by: list[str] = Field(default_factory=list, description="List of fields to supergroup the bigwigs")

    @field_validator("directory")
    def validate_directory(cls, v: str) -> str:
        """
        Make sure that at least the parent directory exists.
        """
        parent = Path(v).parent
        if not parent.is_dir():
            raise ValueError(f"Parent directory {parent} does not exist.")
        return v

    @field_validator("name")
    def validate_name(cls, v: str) -> str:
        if not v:
            raise ValueError("Name must not be empty.")

        if not re.match(r"^[a-zA-Z0-9_.-]+$", v):
            raise ValueError("Name must only contain alphanumeric characters, underscores, dashes, and periods.")
        
        return v

    @classmethod
    def for_assay(cls, assay: Assay) -> "UCSCHubConfig":
        
        match assay:
            case Assay.RNA:
                return cls(
                    subgroup_by=["method", "norm", "strand"],
                    overlay_by=["samplename", "method", "norm"],
                )
            
            case Assay.MCC:
                return cls(
                    color_by = ['viewpoint'],
                    subgroup_by = ['norm', 'viewpoint']
                )
            
            case _:
                return cls()    

class RNAQuantificationConfig(BaseModel, PathValidatorMixin):
    """Configuration for RNA quantification."""

    method: QuantificationMethod
    salmon_index: str | None = None
    run_deseq2: bool = False

    @field_validator("salmon_index")
    def validate_salmon_index(cls, v: str | None) -> str | None:
        return cls.validate_path_exists(v, "Salmon index path")


class SNPCallingConfig(BaseModel, PathValidatorMixin):
    """Configuration for SNP calling."""

    method: SNPCallingMethod
    snp_database: str | None = None

    @field_validator("snp_database")
    def validate_snp_database(cls, v: str | None, info) -> str | None:
        return cls.validate_required_when(v, False, "snp_database")
    
    @computed_field
    @property
    def annotate_snps(self) -> bool:
        """Whether to annotate SNPs (computed from snp_database presence)."""
        return getattr(self, "snp_database", None) is not None


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

    method: MethylationMethod



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

