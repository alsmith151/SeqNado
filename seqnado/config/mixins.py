
from pathlib import Path
from typing import Union, Optional, Callable, Any
from enum import Enum
from pydantic import BaseModel, computed_field, field_validator, ValidationInfo, model_validator


class ComputedFieldOverrideMixin(BaseModel):
    """Mixin class that allows computed fields to be overridden via explicit values.
    
    When a config is loaded from YAML with explicit values for computed fields,
    those values take precedence over the computed logic.
    
    Usage:
        1. Inherit from this mixin in your pydantic model
        2. Define __computed_fields__ as a set of field names (class attribute)
        3. Use _get_or_compute(field_name, compute_fn) in your computed_field properties
    """
    
    @model_validator(mode='before')
    @classmethod
    def _capture_computed_field_overrides(cls, values: dict[str, Any]) -> dict[str, Any]:
        """Capture explicit values for computed fields before validation."""
        if not isinstance(values, dict):
            return values
            
        # Get computed field names from class attribute
        computed_field_names = getattr(cls, '__computed_fields__', set())
        if not computed_field_names:
            return values
        
        # Store explicit overrides in a private dict
        explicit_overrides = {}
        for field_name in computed_field_names:
            if field_name in values:
                explicit_overrides[field_name] = values.pop(field_name)
        
        if explicit_overrides:
            values['_explicit_overrides'] = explicit_overrides
            
        return values
    
    def _get_or_compute(self, field_name: str, compute_fn: Callable[[], Any]) -> Any:
        """Get explicit override value if present, otherwise compute the value.
        
        Args:
            field_name: Name of the computed field
            compute_fn: Function to compute the value if no override exists
            
        Returns:
            The explicit override value or the computed value
        """
        explicit_overrides = getattr(self, '_explicit_overrides', {})
        if field_name in explicit_overrides:
            return explicit_overrides[field_name]
        return compute_fn()
    
    class Config:
        # Allow the _explicit_overrides private field
        extra = 'allow'


class CommonComputedFieldsMixin(ComputedFieldOverrideMixin):
    """Mixin class providing common computed fields for assay configurations."""
    
    __computed_fields__ = {
        'create_bigwigs', 
        'plot_with_plotnado', 
        'create_dataset', 
        'create_ucsc_hub', 
        'has_spikein'
    }

    @computed_field
    @property
    def create_bigwigs(self) -> bool:
        """Whether to make bigwigs 
        (respects explicit config, else computed from bigwigs config presence).
        Returns False for methylation and SNP assays."""
        def compute():
            # Prevent bigwigs for meth and snp assays
            if isinstance(self, (MethylationMixin, SNPCallingMixin)):
                return False
            # Otherwise, fallback to computed logic
            return getattr(self, "bigwigs", None) is not None
        
        return self._get_or_compute('create_bigwigs', compute)

    @computed_field
    @property
    def plot_with_plotnado(self) -> bool:
        """Whether to perform plotting (computed from plotting config presence)."""
        return self._get_or_compute(
            'plot_with_plotnado',
            lambda: getattr(self, "plotting", None) is not None
        )

    @computed_field
    @property
    def create_dataset(self) -> bool:
        """Whether to make dataset."""
        return self._get_or_compute(
            'create_dataset',
            lambda: getattr(self, "dataset_for_ml", None) is not None
        )

    @computed_field
    @property
    def create_ucsc_hub(self) -> bool:
        """Whether to make UCSC hub (computed from ucsc_hub config presence)."""
        return self._get_or_compute(
            'create_ucsc_hub',
            lambda: getattr(self, "ucsc_hub", None) is not None
        )
    
    @computed_field
    def has_spikein(self) -> bool:
        """Whether to use spike-in normalization (computed from spikein config presence)."""
        return self._get_or_compute(
            'has_spikein',
            lambda: getattr(self, "spikein", None) is not None
        )
    
    @model_validator(mode='after')
    def _validate_computed_field_consistency(self) -> 'CommonComputedFieldsMixin':
        """Validate that explicit computed field overrides are consistent with config."""
        explicit_overrides = getattr(self, '_explicit_overrides', {})
        
        # Validate create_bigwigs
        if explicit_overrides.get('create_bigwigs') is True:
            if getattr(self, 'bigwigs', None) is None:
                raise ValueError(
                    "create_bigwigs=True requires 'bigwigs' configuration to be defined"
                )
        
        # Validate plot_with_plotnado
        if explicit_overrides.get('plot_with_plotnado') is True:
            if getattr(self, 'plotting', None) is None:
                raise ValueError(
                    "plot_with_plotnado=True requires 'plotting' configuration to be defined"
                )
        
        # Validate create_dataset
        if explicit_overrides.get('create_dataset') is True:
            if getattr(self, 'dataset_for_ml', None) is None:
                raise ValueError(
                    "create_dataset=True requires 'dataset_for_ml' configuration to be defined"
                )
        
        # Validate create_ucsc_hub
        if explicit_overrides.get('create_ucsc_hub') is True:
            if getattr(self, 'ucsc_hub', None) is None:
                raise ValueError(
                    "create_ucsc_hub=True requires 'ucsc_hub' configuration to be defined"
                )
        
        # Validate has_spikein
        if explicit_overrides.get('has_spikein') is True:
            if getattr(self, 'spikein', None) is None:
                raise ValueError(
                    "has_spikein=True requires 'spikein' configuration to be defined"
                )
        
        return self
    
    @field_validator("create_geo_submission_files", mode="before", check_fields=False)
    def validate_geo_submission_files(cls, v):
        """Ensure geo submission files are created only if the config is set."""
        return v if v else False


class PeakCallingMixin(ComputedFieldOverrideMixin):
    """Mixin for assays that support peak calling."""
    
    __computed_fields__ = {'call_peaks'}

    @computed_field
    @property
    def call_peaks(self) -> bool:
        """Whether to call peaks (computed from peak_calling config presence)."""
        return self._get_or_compute(
            'call_peaks',
            lambda: getattr(self, "peak_calling", None) is not None
        )
    
    @model_validator(mode='after')
    def _validate_peak_calling_consistency(self) -> 'PeakCallingMixin':
        """Validate that explicit call_peaks override is consistent with config."""
        explicit_overrides = getattr(self, '_explicit_overrides', {})
        
        if explicit_overrides.get('call_peaks') is True:
            if getattr(self, 'peak_calling', None) is None:
                raise ValueError(
                    "call_peaks=True requires 'peak_calling' configuration to be defined"
                )
        
        return self


class SNPCallingMixin(ComputedFieldOverrideMixin):
    """Mixin for assays that support SNP calling."""
    
    __computed_fields__ = {'call_snps'}

    @computed_field
    @property
    def call_snps(self) -> bool:
        """Whether to call SNPs (computed from snp_calling config presence)."""
        return self._get_or_compute(
            'call_snps',
            lambda: getattr(self, "snp_calling", None) is not None
        )
    
    @model_validator(mode='after')
    def _validate_snp_calling_consistency(self) -> 'SNPCallingMixin':
        """Validate that explicit call_snps override is consistent with config."""
        explicit_overrides = getattr(self, '_explicit_overrides', {})
        
        if explicit_overrides.get('call_snps') is True:
            if getattr(self, 'snp_calling', None) is None:
                raise ValueError(
                    "call_snps=True requires 'snp_calling' configuration to be defined"
                )
        
        return self


class MethylationMixin(ComputedFieldOverrideMixin):
    """Mixin for assays that support methylation calling."""
    
    __computed_fields__ = {'call_methylation'}

    @computed_field
    @property
    def call_methylation(self) -> bool:
        """Whether to call methylation (computed from methylation config presence)."""
        return self._get_or_compute(
            'call_methylation',
            lambda: getattr(self, "methylation", None) is not None
        )
    
    @model_validator(mode='after')
    def _validate_methylation_consistency(self) -> 'MethylationMixin':
        """Validate that explicit call_methylation override is consistent with config."""
        explicit_overrides = getattr(self, '_explicit_overrides', {})
        
        if explicit_overrides.get('call_methylation') is True:
            if getattr(self, 'methylation', None) is None:
                raise ValueError(
                    "call_methylation=True requires 'methylation' configuration to be defined"
                )
        
        return self


class PathValidatorMixin:
    """Mixin providing common path validation methods."""

    @staticmethod
    def validate_path_exists(
        v: Path | str | None, field_name: str = "path", info: ValidationInfo | None = None
    ) -> Path | str | None:
        """Validate that a path exists if provided.

        Skips validation if 'skip_path_validation' is True in the validation context.
        This allows creating template configs with placeholder paths.
        """
        if v is not None:
            # Skip validation if context flag is set (for template generation)
            if info and info.context and info.context.get("skip_path_validation", False):
                return v
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