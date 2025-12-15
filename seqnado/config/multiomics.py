from pathlib import Path
from typing import Annotated

from pydantic import BaseModel, BeforeValidator, Field


def none_str_to_none(v):
    """Convert string 'none' to None."""
    if isinstance(v, str) and v.strip().lower() == "none":
        return None
    return v


class MultiomicsConfig(BaseModel):
    """Configuration for multiomics analysis combining multiple assays.

    This config can be used in two modes:
    1. File-based mode: Only stores multiomics settings, assay configs discovered from files
    2. Unified mode: Stores all assay configs and multiomics settings in one place
    """

    # List of assay names to include (used for validation and config generation)
    assays: list[str] = Field(
        default_factory=list,
        description="List of assay names to include in multiomics analysis (e.g., ['atac', 'chip', 'rna'])",
    )

    # Optional: store full assay configs (for unified config mode)
    # This allows a single config_multiomics.yaml to contain everything
    assay_configs: dict[str, dict] = Field(
        default_factory=dict,
        description="Dictionary mapping assay names to their SeqnadoConfig dicts (optional)",
    )

    output_dir: str = Field(
        default="seqnado_output/", description="Output directory for multiomics results"
    )

    # Multiomics-specific analysis options
    create_heatmaps: bool = Field(
        default=True, description="Generate heatmaps for multiomics data"
    )
    create_dataset: bool = Field(
        default=True, description="Generate ML-ready dataset combining all assays"
    )
    create_summary: bool = Field(
        default=True, description="Generate summary report of multiomics analysis"
    )

    # Optional settings for multiomics analysis
    regions_bed: Annotated[Path | None, BeforeValidator(none_str_to_none)] = Field(
        default=None,
        description="BED file with regions of interest for multiomics analysis",
    )
    binsize: int | None = Field(
        default=None, description="Bin size for genome-wide multiomics analysis"
    )
