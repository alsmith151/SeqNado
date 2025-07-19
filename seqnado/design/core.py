"""Core enums, constants, and utility functions for SeqNado design module."""

import pathlib
import re
from enum import Enum
from typing import Optional, Union
from pydantic import BaseModel, Field, computed_field, field_validator


# =============================================================================
# ENUMS
# =============================================================================

class Assay(Enum):
    """Supported sequencing assay types."""
    RNA = "rna"
    ATAC = "atac"
    SNP = "snp"
    CHIP = "chip"
    CAT = "cat"
    METH = "meth"
    MCC = "mcc"
    CRISPR = "crispr"
    
    @classmethod
    def non_ip_assays(cls):
        """Return assays that don't require IP (immunoprecipitation)."""
        ip_assays = {cls.CHIP, cls.CAT}
        return [assay for assay in cls if assay not in ip_assays]
    
    @classmethod
    def ip_assays(cls):
        """Return assays that require IP (immunoprecipitation)."""
        return [cls.CHIP, cls.CAT]


class PileupMethod(Enum):
    """Methods for creating pileup files."""
    DEEPTOOLS = "deeptools"
    HOMER = "homer"
    BAMNADO = "bamnado"


class ScaleMethod(Enum):
    """Methods for scaling genomic data."""
    UNSCALED = "unscaled"
    CSAW = "csaw"
    CPM = "cpm"
    RPKM = "rpkm"
    SPIKEIN = "spikein"
    MERGED = "merged"


class PeakCallingMethod(Enum):
    """Methods for calling peaks."""
    MACS2 = "macs2"
    MACS3 = "macs3"
    HOMER = "homer"
    LANCEOTRON = "lanceotron"
    SEACR = "seacr"
    LANCEOTRON_MCC = "lanceotron-mcc"


# =============================================================================
# CONSTANTS
# =============================================================================
ILLUMINA_FILENAME_PATTERNS = {
    r"_S\d+_": "_",
    r"_L00\d_": "_",
    r"_R?[12](_001)?$": "_",
    r"__": "_",
    r"_$": "",
}

INPUT_CONTROL_SUBSTRINGS = ["input", "mock", "igg", "control"]

NONE_VALUES = [
    None, "None", "none", "null", "Null", "NULL", 
    ".", "", "NA"
]




# =============================================================================
# Models
# =============================================================================

class Metadata(BaseModel):
    """Metadata for samples. Optional fields can be set to None."""
    assay: Assay | None = Field(
        default=None,
        description="Assay type, should be one of the Assay enum values"
    )
    consensus_group: str | None = Field(
        default=None,
        description="Grouping variable for merging samples into consensus, can be None if not applicable"
    )
    norm_group: str = Field(
        default="all",
        description="Grouping variable for scaling samples, defaults to 'all' if not specified"
    )
    deseq2: str | None = Field(
        default=None,
        description="DESeq2 metadata for sample, can be None if not applicable"
    )

    @field_validator("deseq2", "consensus_group")
    @classmethod
    def prevent_none(cls, v):
        import numpy as np
        import pandas as pd
        """Ensure metadata fields are not set to None."""
        # Define a list of values that should be treated as 'None'
        # This includes None, pd.NA, np.nan, and common string representations of None

        none_vals = [
            None,
            "None",
            "none",
            "null",
            "Null",
            "NULL",
            ".",
            "",
            "NA",
            np.nan,
        ]
        if any([v == n for n in none_vals]):
            assert v is not None, "None is not allowed when setting metadata"

        # Check if it a pandas NA
        if v is pd.NA:
            raise ValueError("None is not allowed when setting metadata")
        return v



# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def predict_organism(genome: str) -> str:
    """Predict organism from genome identifier."""
    if "hg" in genome:
        return "Homo sapiens"
    elif "mm" in genome:
        return "Mus musculus"
    else:
        return "Unknown"


def is_valid_path(path: Optional[Union[str, pathlib.Path]]) -> bool:
    """Check if a path is valid and exists."""
    if path is None:
        return False
    
    try:
        p = pathlib.Path(path)
        return p.exists() and str(p) not in ["-", ".", "", "None"]
    except (TypeError, ValueError):
        return False


# =============================================================================
# STRING/FILENAME PROCESSING FUNCTIONS
# =============================================================================

def clean_sample_name(name: str, patterns: dict = ILLUMINA_FILENAME_PATTERNS) -> str:
    """Clean sample name using regex patterns."""
    for pattern, replacement in patterns.items():
        name = re.sub(pattern, replacement, name)
    return name


def extract_read_number(filename: str) -> Optional[int]:
    """Extract read number from FASTQ filename."""
    regex_std_illumina_paired = re.compile(r".*_R?([12])(_001)?")
    
    match = regex_std_illumina_paired.match(filename)
    if match:
        return int(match.group(1))
    return None


def is_control_sample(ip_name: str) -> bool:
    """Check if sample is a control based on IP name."""
    return any(substring in ip_name.lower() for substring in INPUT_CONTROL_SUBSTRINGS)
