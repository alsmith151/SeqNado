import pandera.pandas as pa
from pandera.typing.pandas import Series
from typing import Optional, Union, Annotated
import pandas as pd
from .core import Assay

AssayCategory = pd.CategoricalDtype(categories=[a.value for a in Assay])

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


class DesignDataFrame(pa.DataFrameModel):
    """Base class for design dataframes with common sample identification."""
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
    r1: Series[str] = pa.Field(coerce=True)
    r2: Series[str] | None = pa.Field(coerce=True, nullable=True, description="Optional second read for paired-end data")
    r1_control: Series[str] | None = pa.Field(
        coerce=True, nullable=True, description="Optional control read 1 for paired-end data"
    )
    r2_control: Series[str] | None = pa.Field(
        coerce=True, nullable=True, description="Optional control read 2 for paired-end data"
    )
    ip: Series[str] | None = pa.Field(coerce=True, nullable=True, description="Optional IP read for IP-seq data")
    control: Series[str] | None = pa.Field(coerce=True, nullable=True, description="Optional control sample name for IP-seq data")
    assay: Series[Annotated[pd.CategoricalDtype, [a.value for a in Assay], True]] | None = pa.Field(
        default=None,
        description="Assay type, should be one of the Assay enum values",
        coerce=True,
        nullable=True,
    )

    @pa.dataframe_check
    def check_sample_name(cls, df: pd.DataFrame) -> Series[bool]:
        """Ensure that either the sample_name or sample_name + 'ip' is unique."""
        
        if 'ip' in df.columns:
            # If 'ip' column exists, check uniqueness of sample_name + ip
            unique_combination = df['sample_name'] + df['ip'].fillna('')
            return unique_combination.is_unique
        else:
            # If 'ip' column does not exist, check uniqueness of sample_name
            return df['sample_name'].is_unique


        

    