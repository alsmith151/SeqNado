from enum import Enum
from typing import Any, List, Literal, Optional, Union

import pandas as pd
from pydantic import BaseModel, Field, computed_field, field_validator

from seqnado import Assay, Molecule, Organism


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
    single_or_paired: Literal["single", "paired-end"] = "paired-end"
    instrument_model: str | None = "NovaSeq X"
    description: str | None = None
    processed_data_file: list[str] = Field(default_factory=list)
    raw_file: list[str] = Field(default_factory=list)

    @computed_field
    @property
    def library_strategy(self) -> str:
        """Return the library strategy based on the assay."""
        return self.assay.value

    @field_validator("antibody", mode="before")
    def validate_antibody_with_assay(cls, v: Optional[str], values: dict[str, Any]) -> Optional[str]:
        """Ensure antibody is set for ChIP and CUT&TAG assays."""
        assay = values.get("assay")
        if assay in Assay.ip_assays() and not v:
            raise ValueError(f"Antibody must be specified for {assay.clean_name} assays.")
        elif assay not in Assay.ip_assays() and v:
            raise ValueError(f"Antibody should not be specified for {assay.clean_name} assays.")
        return v

    @computed_field
    @property
    def molecule(self) -> Molecule:
        """Return the molecule type based on the assay."""
        match self.assay:
            case Assay.RNA:
                if any(n in self.title.lower() for n in ["tt-seq", "nasc", "point"]):
                    return Molecule.rna_nuclear
                return Molecule.rna_polya
            case _:
                return Molecule.dna_genomic

    def to_series(self) -> pd.Series:
        """Convert the GEOSample to a pandas Series with dynamic file columns."""

        # Static fields
        data = self.model_dump(exclude_none=True)

        if self.assay in Assay.ip_assays():
            data['ChIP antibody'] = self.antibody
            del data['antibody']
    
        # Add processed data files dynamically
        processed_files = self.processed_data_file
        for i, file in enumerate(processed_files):
            key = "processed data file" if i == 0 else f"processed data file {i}"
            data[key] = file

        # Add raw files dynamically
        raw_files = self.raw_file
        for i, file in enumerate(raw_files):
            key = "raw file" if i == 0 else f"raw file {i}"
            data[key] = file

        return pd.Series(data)
    
    @classmethod
    def from_series(cls, series: pd.Series) -> "GEOSample":
        """Create a GEOSample instance from a pandas Series."""
        data = series.to_dict()

        # Extract dynamic file columns
        processed_files = [data.pop(key) for key in list(data.keys()) if key.startswith("processed data file")]
        raw_files = [data.pop(key) for key in list(data.keys()) if key.startswith("raw file")]

        # Create the GEOSample instance
        return cls(
            **data,
            processed_data_file=processed_files,
            raw_file=raw_files
        )


class GEOSamples(BaseModel):
    samples: list[GEOSample]

    def to_dataframe(self):
        df = pd.concat([s.to_series() for s in self.samples], axis=1).T
        df.columns = df.columns.str.replace(r"\s\d+$", "", regex=True).str.strip()
        return df 
   


        
        
