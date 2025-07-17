"""FASTQ file handling classes for genomic sequencing workflows."""

from __future__ import annotations

import pathlib
import re

import numpy as np
import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field, computed_field, field_validator

from .core import (
    clean_sample_name,
    extract_read_number,
    is_control_sample,
    is_valid_path,
)


# =============================================================================
# BASE FASTQ CLASSES
# =============================================================================

class FastqFile(BaseModel):
    """Represents a single FASTQ file with metadata extraction."""
    
    path: pathlib.Path
    use_resolved_name: bool = False

    def model_post_init(self, __context: dict[str, any] | None) -> None:
        """Validate file path after initialization."""
        if self.use_resolved_name:
            self.path = self.path.resolve()
        else:
            self.path = self.path.absolute()

        if not is_valid_path(self.path):
            raise FileNotFoundError(f"FASTQ file not found: {self.path}")

    @property
    def stem(self) -> str:
        """Get file stem without .gz extension."""
        return pathlib.Path(str(self.path).removesuffix(".gz")).stem

    @computed_field
    @property
    def sample_name(self) -> str:
        """Extract clean sample name from filename."""
        name = self.stem
        if name.endswith("_001"):
            name = name.removesuffix("_001")
        return name

    @computed_field
    @property
    def sample_base(self) -> str:
        """Extract base sample name by cleaning Illumina patterns."""
        return clean_sample_name(self.sample_name)

    @computed_field
    @property
    def read_number(self) -> int | None:
        """Extract read number (1 or 2) from filename."""
        read_num = extract_read_number(self.sample_name)
        if read_num is None:
            logger.debug(
                f"Could not extract read number from '{self.sample_name}', "
                f"assuming single-end sequencing"
            )
        return read_num

    @computed_field
    @property
    def is_paired(self) -> bool:
        """Check if file is part of a paired-end experiment."""
        return self.read_number is not None

    @computed_field
    @property
    def is_lane(self) -> bool:
        """Check if file is lane-split."""
        return "_L00" in self.sample_name

    def __lt__(self, other: FastqFile) -> bool:
        """Compare FastqFile objects by path."""
        return self.path < other.path

    def __gt__(self, other: FastqFile) -> bool:
        """Compare FastqFile objects by path."""
        return self.path > other.path

    def __eq__(self, other: object) -> bool:
        """Check equality of FastqFile objects by path."""
        if not isinstance(other, FastqFile):
            return NotImplemented
        return self.path == other.path

    def __hash__(self) -> int:
        """Make FastqFile hashable based on path."""
        return hash(self.path)


# =============================================================================
# IP EXPERIMENT CLASSES
# =============================================================================

class FastqFileIP(FastqFile):
    """FASTQ file for IP experiments with antibody information."""
    
    ip: str = Field(default=None, description="IP antibody performed on the sample")
    is_control: bool = Field(default=None, description="Whether sample is a control")

    def model_post_init(self, __context: dict[str, any] | None) -> None:
        """Initialize and predict IP information."""
        super().model_post_init(__context)
        
        if self.ip is None:
            self.ip = self._predict_ip()

        if self.is_control is None:
            self.is_control = self._predict_is_control()

    def _predict_ip(self) -> str | None:
        """Predict IP antibody from sample name."""
        try:
            return self.sample_base.split("_")[-1]
        except IndexError:
            logger.warning(f"Could not predict IP for {self.sample_base}")
            return None

    def _predict_is_control(self) -> bool:
        """Predict if sample is a control."""
        if not self.ip:
            return False
        return is_control_sample(self.ip)

    @computed_field
    @property
    def sample_base_without_ip(self) -> str:
        """Get sample base name without antibody suffix."""
        if not self.ip:
            return self.sample_base
            
        # Remove IP and common Illumina suffixes
        pattern = rf"(_{re.escape(self.ip)})?(_S\d+)?(_L00\d)?(_R?[12])?(_001)?"
        return re.sub(pattern, "", self.sample_name)

    @field_validator("ip")
    @classmethod
    def allow_na_or_nan(cls, v: str | None | float) -> str | None | float:
        """Allow None, pd.NA, or np.nan values for IP."""
        if v is None or v is pd.NA or (isinstance(v, float) and np.isnan(v)):
            return v
        if not isinstance(v, str):
            raise ValueError("ip must be a string, None, pd.NA, or np.nan")
        return v


# =============================================================================
# FASTQ SET CLASSES
# =============================================================================

class FastqSet(BaseModel):
    """Represents a set of FASTQ files for a single sample."""
    
    name: str = Field(default=None, description="Sample name")
    r1: FastqFile
    r2: FastqFile | None = None

    @property
    def fastq_paths(self) -> list[pathlib.Path]:
        """Get list of FASTQ file paths."""
        paths = [self.r1.path]
        if self.is_paired and self.r2:
            paths.append(self.r2.path)
        return paths

    @property
    def is_paired(self) -> bool:
        """Check if sample has paired-end reads."""
        return self.r2 is not None and self.r2.path.exists()

    @classmethod
    def from_fastq_files(
        cls, 
        fq: list[FastqFile], 
        **kwargs
    ) -> FastqSet:
        """Create FastqSet from list of FastqFile objects."""
        if not fq:
            raise ValueError("No FASTQ files provided")
            
        sample_name = fq[0].sample_base

        match len(fq):
            case 1:
                return cls(name=sample_name, r1=fq[0], **kwargs)
            case 2:
                return cls(name=sample_name, r1=fq[0], r2=fq[1], **kwargs)
            case _:
                raise ValueError(
                    f"Invalid number of FASTQ files for {sample_name}: {len(fq)}. "
                    f"Expected 1 or 2 files."
                )


class FastqSetIP(FastqSet):
    """FASTQ set for IP experiments."""
    
    r1: FastqFileIP
    r2: FastqFileIP | None = None
    ip: str = Field(default=None, description="IP antibody performed")

    @property
    def ip_or_control_name(self) -> str:
        """Get IP or control name."""
        return self.ip if self.ip else self.r1.ip

    @property
    def sample_name(self) -> str:
        """Get full sample name including IP."""
        return f"{self.name}_{self.ip_or_control_name}"

    @property
    def sample_name_base(self) -> str:
        """Get base sample name without IP."""
        return self.r1.sample_base_without_ip

    @property
    def is_control(self) -> bool:
        """Check if this is a control sample."""
        return self.r1.is_control
