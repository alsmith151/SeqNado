from typing import Any, List, Protocol, Optional
from pathlib import Path
from itertools import chain

from pydantic import BaseModel, computed_field, Field, field_validator
from typing import List, Literal
from snakemake.io import expand
from loguru import logger

from seqnado import (
    Assay,
    PileupMethod,
    ScaleMethod,
    AssaysWithPeakCalling,
    PeakCallingMethod,
    MethylationMethod,
    QuantificationMethod,
)
from seqnado.core import AssaysWithHeatmaps, AssaysWithSpikein
from seqnado.inputs import SampleCollection, IPSampleCollection, SampleGroups
from seqnado.config import SeqnadoConfig


class FileCollection(Protocol):
    @property
    def files(self) -> List[str]:
        """Return a list of file paths."""
        pass


class QCFiles(BaseModel):
    assay: Assay
    samples: SampleCollection | IPSampleCollection

    @property
    def default_files(self) -> list[str]:
        return [
            "seqnado_output/seqnado_report.html",
        ]

    @property
    def qualimap_files(self) -> list[str]:
        match self.assay:
            case Assay.RNA:
                return expand(
                    "seqnado_output/qc/qualimap_rnaseq/{sample}/qualimapReport.html",
                    sample=self.samples.sample_names,
                )
            case _:
                return expand(
                    "seqnado_output/qc/qualimap_bamqc/{sample}/qualimapReport.html",
                    sample=self.samples.sample_names,
                )

    @computed_field
    @property
    def files(self) -> List[str]:
        return [*self.default_files, *self.qualimap_files]


class BigWigFiles(BaseModel):
    assay: Assay
    names: list[str] = Field(default_factory=list)
    pileup_methods: list[PileupMethod]
    scale_methods: list[ScaleMethod] = [ScaleMethod.UNSCALED]
    prefix: Path | None = "seqnado_output/bigwigs/"

    @property
    def is_rna(self) -> bool:
        return self.assay == Assay.RNA

    @property
    def incompatible_methods(self) -> dict[PileupMethod, list[ScaleMethod]]:
        return {
            PileupMethod.HOMER: [ScaleMethod.CSAW, ScaleMethod.SPIKEIN],
            PileupMethod.BAMNADO: [ScaleMethod.CSAW, ScaleMethod.SPIKEIN],
        }

    def _is_compatible(self, method: PileupMethod, scale: ScaleMethod) -> bool:
        return scale not in self.incompatible_methods.get(method, [])

    def generate_bigwig_paths(self) -> list[str]:
        paths = []

        for method in self.pileup_methods:
            for scale in self.scale_methods:
                if not self._is_compatible(method, scale):
                    continue

                if self.is_rna:
                    for name in self.names:
                        for strand in ["plus", "minus"]:
                            path = f"{self.prefix}{method.value}/{scale.value}/{name}_{strand}.bigWig"
                            paths.append(path)
                else:
                    for name in self.names:
                        path = (
                            f"{self.prefix}{method.value}/{scale.value}/{name}.bigWig"
                        )
                        paths.append(path)

        return paths

    @computed_field
    @property
    def files(self) -> list[str]:
        return self.generate_bigwig_paths()


class PeakCallingFiles(BaseModel):
    assay: AssaysWithPeakCalling
    names: list[str]
    peak_calling_method: list[PeakCallingMethod]
    prefix: Optional[str] = "seqnado_output/peaks/"

    @property
    def peak_files(self) -> list[str]:
        return expand(
            self.prefix + "{method}/{sample}.bed",
            sample=self.names,
            method=[m.value for m in self.peak_calling_method],
        )

    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of peak calling files."""
        return self.peak_files


class HeatmapFiles(BaseModel):
    assay: AssaysWithHeatmaps

    @property
    def heatmap_files(self) -> list[str]:
        return [
            "seqnado_output/heatmap/heatmap.pdf",
            "seqnado_output/heatmap/metaplot.pdf",
        ]

    @computed_field
    @property
    def files(self) -> list[str]:
        return self.heatmap_files


class HubFiles(BaseModel):
    hub_dir: Path
    hub_name: str

    @property
    def hub_txt(self) -> Path:
        return self.hub_dir / f"{self.hub_name}.hub.txt"

    @computed_field
    @property
    def files(self) -> List[str]:
        return [str(self.hub_txt)]


class SpikeInFiles(BaseModel):
    assay: AssaysWithSpikein
    names: list[str]

    @property
    def norm_factors(self):
        return "seqnado_output/resources/normalisation_factors.tsv"

    @computed_field
    @property
    def files(self) -> list[str]:
        return [self.norm_factors]


class PlotFiles(BaseModel):
    coordinates: Path
    file_format: Literal["svg", "png", "pdf"] = "svg"

    @property
    def plot_names(self):
        import pandas as pd

        plots = []

        try:
            # Read BED file using pandas (pyranges replacement)
            bed_columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
            coords_df = pd.read_csv(
                str(self.coordinates), sep="\t", header=None, comment="#"
            )
            coords_df.columns = bed_columns[: len(coords_df.columns)]

            outdir = Path("seqnado_output/genome_browser_plots/")
            for region in coords_df.itertuples():
                fig_name = (
                    f"{region.Chromosome}-{region.Start}-{region.End}"
                    if not hasattr(region, "Name") or not region.Name
                    else region.Name
                )
                plots.append(outdir / f"{fig_name}.{self.file_format}")

        except FileNotFoundError:
            logger.warning(
                f"Could not find plotting coordinates file: {self.coordinates}"
            )

        return plots

    @computed_field
    @property
    def files(self) -> List[str]:
        """Return a list of plot files."""
        return self.plot_names

class SNPFilesRaw(BaseModel):
    assay: Assay
    names: list[str]

    @property
    def snp_files(self) -> list[str]:
        return expand(
            "seqnado_output/variant/{sample}.vcf.gz",
            sample=self.names,
        )
    
    @computed_field
    @property
    def files(self) -> list[str]:
        return self.snp_files

class SNPFilesAnnotated(BaseModel):
    assay: Assay
    names: list[str]

    @property
    def anno_snp_files(self) -> list[str]:
        return expand(
            "seqnado_output/variant/{sample}.anno.vcf.gz",
            sample=self.names,
        )

    @computed_field
    @property
    def files(self) -> list[str]:
        return self.anno_snp_files


class MethylationFiles(BaseModel):
    assay: Assay
    names: list[str]
    genomes: List[str]
    method: MethylationMethod
    prefix: Path | None = "seqnado_output/methylation/"

    @property
    def split_bams_files(self) -> List[str]:
        return expand(
            "seqnado_output/aligned/spikein/{sample}_{genome}.bam",
            sample=self.names,
            genome=self.genomes,
        )
    
    @property
    def methyldackel_files(self) -> List[str]:

        file_pattern = "seqnado_output/methylation/methyldackel/{sample}_{genome}_|METHOD|CpG.bedGraph"
        file_pattern = file_pattern.replace("|METHOD|", self.method.value)
        return expand(
            file_pattern,
            sample=self.names,
            genome=self.genomes,
        )
    
    @property
    def methylation_bias(self) -> List[str]:
        """Return the methylation bias files."""
        files = []
        files.append("seqnado_output/methylation/methylation_conversion.tsv")
        files.extend(
            expand(
                "seqnado_output/methylation/methyldackel/bias/{sample}_{genome}.txt",
                sample=self.names,
                genome=self.genomes,
            )
        )
        return files

    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of methylation files."""
        return [*self.split_bams_files, *self.methyldackel_files, *self.methylation_bias]
    


class BigBedFiles(BaseModel):
    bed_files: list[Path] = Field(default_factory=list)

    @field_validator("bed_files", mode="before")
    def validate_bed_files(cls, v: list[Path | str]) -> list[Path]:
        return [Path(f) for f in v]
    
    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of bigBed files."""
        return [str(f.with_suffix(".bb")) for f in self.bed_files if f.suffix == ".bed"]


class ContactFiles(BaseModel):
    assay: Assay
    names: list[str]
    prefix: Path | None = "seqnado_output/contacts/"

    @property
    def cooler_files(self) -> List[str]:
        return expand(
            str(self.prefix).rstrip("/") + "/{group}/{group}.mcool",
            group=self.design_dataframe["merge"].unique().tolist(),
        )

    @property
    def pairs(self) -> List[str]:
        return expand(
            str(self.prefix).rstrip("/") + "/{group}/ligation_junctions/{viewpoint}.pairs.gz",
            group=self.design_dataframe["merge"].unique().tolist(),
            viewpoint=self.viewpoints_grouped,
        )

    @computed_field
    @property
    def files(self) -> List[str]:
        """Return a list of contact files."""
        return [*self.cooler_files, *self.pairs]
    


class QuantificationFiles(BaseModel):
    """Base class for quantification files."""
    
    assay: Assay
    methods: list[QuantificationMethod] = Field(default_factory=list)
    names: list[str]
    groups: SampleGroups
    prefix: Path | None = "seqnado_output/quantification"

    @field_validator("methods", mode="before")
    def validate_methods_and_assays(cls, v: list[QuantificationMethod], values: dict[str, Any]) -> list[QuantificationMethod]:
        assay = values.get("assay")
        if assay == Assay.RNA:
            return [m for m in v if m in [QuantificationMethod.FEATURE_COUNTS, QuantificationMethod.SALMON]]
        return [m for m in v if m == QuantificationMethod.FEATURE_COUNTS]

    @property
    def combined_counts_file(self) -> list[str]:
        """Return the combined read counts file."""
        return expand(self.prefix + "/{methods}/read_counts.tsv", methods=self.methods)
    
    @property
    def grouped_counts_files(self) -> list[str]:
        """Return the grouped read counts files."""
        files = []
        for group in self.groups.groups:
            files.extend(
                expand(
                    f"{self.prefix}{group}/read_counts.tsv",
                    group=group,
                    methods=self.methods,
                )
            )
        return files
    
    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of quantification files."""
        return [*self.combined_counts_file, *self.grouped_counts_files]
            
    

    





    
    
    