from enum import Enum
from typing import List, Optional, Union, Literal
from pydantic import BaseModel, Field, field_validator, computed_field

from seqnado import Assay, Organism


class GEOSample(BaseModel):
    assay: Assay
    library_name: str
    title: str
    organism: Organism
    cell_line: str | None = None
    cell_type: str | None = None
    antibody: str | None = None
    genotype: str | None = None
    treatment: str | None = None
    time: str | None = None
    single_or_paired: Literal["single", "paired-end"]
    instrument_model: str
    description: str | None = None
    processed_data_file: list[str] = Field(default_factory=list)
    raw_file: list[str] = Field(default_factory=list)