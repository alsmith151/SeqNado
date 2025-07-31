
from pathlib import Path
from typing import Union, Optional
from enum import Enum
from pydantic import BaseModel, computed_field, field_validator


class CommonComputedFieldsMixin:
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
    
    @computed_field
    def has_spikein(self) -> bool:
        """Whether to use spike-in normalization (computed from spikein config presence)."""
        return getattr(self, "spikein", None) is not None
    
    @field_validator("create_geo_submission_files", mode="before")
    def validate_geo_submission_files(cls, v):
        """Ensure geo submission files are created only if the config is set."""
        return v if v else False


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