from pathlib import Path
from typing import Annotated

from pydantic import BaseModel, BeforeValidator, Field

from seqnado import Assay


def none_str_to_none(v):
    """Convert string 'none' to None."""
    if isinstance(v, str) and v.strip().lower() == "none":
        return None
    return v


class MultiomicsConfig(BaseModel):
    """Configuration for multiomics analysis combining multiple assays."""

    assays: list[Assay] = Field(
        default_factory=list,
        description="List of assays to include in multiomics analysis",
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
