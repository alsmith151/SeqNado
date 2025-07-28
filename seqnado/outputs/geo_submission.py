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
        """Convert the GEOSample to a pandas Series."""
        

        model_data = self.model_dump(exclude_none=True)        
        processed_data = {
            f"processed data file {i}": f
            for i, f in enumerate(self.processed_data_file)
        }
        raw_data = {f"raw file {i}": f for i, f in enumerate(self.raw_file)}

        if self.assay in Assay.ip_assays():
            model_data["ChIP antibody"] = self.antibody
            del model_data["antibody"]
        
        return pd.Series(
            {
                **model_data,
                **processed_data,
                **raw_data,
                "molecule": self.molecule.value,
            }
        )
