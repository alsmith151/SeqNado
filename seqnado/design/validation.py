import pandera.pandas as pa
from pandera.typing.pandas import Series
from typing import Optional


class ViewpointsFile(pa.DataFrameModel):
    Chromosome: Series[str] = pa.Field(coerce=True)
    Start: Series[int] = pa.Field(coerce=True)
    End: Series[int] = pa.Field(coerce=True)
    Name: Series[str] = pa.Field(coerce=True)
    Strand: Series[str] | None = pa.Field(coerce=True)
    Score: Series[float] | None = pa.Field(coerce=True, nullable=True)

    # Validate the viewpoint names column
    @pa.check("Name")
    def check_viewpoint_names(cls, s: Series[str]) -> Series[bool]:
        # Check that the names do not contain spaces or special characters
        allowed_chars = r"^[a-zA-Z0-9_]+$"

        return s.str.match(allowed_chars)



class DataFrameDesignBase(pa.DataFrameModel):
    """Base class for all design dataframes with common sample identification."""
    sample_name: Series[str] = pa.Field(coerce=True)
    scale_group: Series[str] | None = pa.Field(coerce=True, default="all", description="Grouping variable for scaling samples")
    merge: Series[str] | None = pa.Field(
        default=None,
        description="Grouping variable for merging samples",
        nullable=False,   
    )
    deseq2: Series[str] | None = pa.Field(
        default=None,
        description="DESeq2 metadata for sample, can be None if not applicable",
        nullable=True,
    )


class DataFrameDesign(DataFrameDesignBase):
    """Standard design for paired-end sequencing data."""
    r1: Series[str] = pa.Field(coerce=True)
    r2: Series[str] = pa.Field(coerce=True, nullable=True)


class DataFrameDesignIP(DataFrameDesignBase):
    """Design for IP-seq experiments with IP and optional control samples."""
    ip: Series[str] = pa.Field(coerce=True)
    control: Optional[Series[str]] = pa.Field(coerce=True, nullable=True)
    ip_r1: Series[str] = pa.Field(coerce=True)
    ip_r2: Series[str] = pa.Field(coerce=True, nullable=True)
    control_r1: Series[str] | None = pa.Field(coerce=True, nullable=True)
    control_r2: Series[str] | None = pa.Field(coerce=True, nullable=True)
    