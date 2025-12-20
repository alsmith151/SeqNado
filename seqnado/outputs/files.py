from typing import Protocol, List
from pathlib import Path

from pydantic import BaseModel, computed_field, Field, field_validator
from typing import Literal
from snakemake.io import expand
from loguru import logger

from seqnado import (
    Assay,
    PileupMethod,
    DataScalingTechnique,
    PeakCallingMethod,
    MethylationMethod,
    QuantificationMethod,
)
from seqnado.core import AssaysWithHeatmaps, AssaysWithSpikein, AssaysWithPeakCalling
from seqnado.inputs import FastqCollection, FastqCollectionForIP, SampleGroups, BamCollection, BigWigCollection
from seqnado.config.configs import QCConfig


class FileCollection(Protocol):
    @property
    def files(self) -> List[str]:
        """Return a list of file paths."""
        pass

class BasicFileCollection(BaseModel):
    files: List[str] = Field(default_factory=list)


class QCFiles(BaseModel):
    assay: Assay
    samples: FastqCollection | FastqCollectionForIP | BamCollection | BigWigCollection
    config: QCConfig = QCConfig()
    output_dir: str = "seqnado_output"

    @property
    def fastqc_files(self) -> list[str]:
        if isinstance(self.samples, (FastqCollection, FastqCollectionForIP)):
            return expand(
                f"{self.output_dir}/qc/fastqc_raw/{{sample}}_{{read}}_fastqc.html",
                sample=self.samples.sample_names,
                read=["1", "2"],
            )
        return []

    @property
    def fastqscreen_files(self) -> list[str]:
        if isinstance(self.samples, (FastqCollection, FastqCollectionForIP)) and self.config.run_fastq_screen:
            return expand(
                f"{self.output_dir}/qc/fastq_screen/{{sample}}_{{read}}_screen.html",
                sample=self.samples.sample_names,
                read=["1", "2"],
            )
        return []
    
    @property
    def qualimap_files(self) -> list[str]:
        if not isinstance(self.samples, (BigWigCollection)):
            match self.assay:
                case Assay.RNA:
                    return expand(
                        f"{self.output_dir}/qc/qualimap_rnaseq/{{sample}}/qualimapReport.html",
                        sample=self.samples.sample_names,
                    )
                case _:
                    return expand(
                        f"{self.output_dir}/qc/qualimap_bamqc/{{sample}}/qualimapReport.html",
                        sample=self.samples.sample_names,
                    )
        return []

    @computed_field
    @property
    def files(self) -> List[str]:
        return [*self.fastqc_files, *self.fastqscreen_files, *self.qualimap_files]
    
class SeqNadoReportFile(BaseModel):
    output_dir: str = "seqnado_output"

    @property
    def report_file(self) -> list[str]:
        return [
            f"{self.output_dir}/seqnado_report.html",
        ]
    
    @computed_field
    @property
    def files(self) -> list[str]:
        return self.report_file


class BigWigFiles(BaseModel):
    assay: Assay
    names: list[str] = Field(default_factory=list)
    pileup_methods: list[PileupMethod]
    scale_methods: list[DataScalingTechnique] = [DataScalingTechnique.UNSCALED]
    output_dir: str = "seqnado_output"

    @property
    def prefix(self) -> str:
        return f"{self.output_dir}/bigwigs/"

    @property
    def is_rna(self) -> bool:
        return self.assay == Assay.RNA

    @property
    def incompatible_methods(self) -> dict[PileupMethod, list[DataScalingTechnique]]:
        return {
            PileupMethod.HOMER: [DataScalingTechnique.CSAW, DataScalingTechnique.SPIKEIN],
            PileupMethod.BAMNADO: [DataScalingTechnique.CSAW, DataScalingTechnique.SPIKEIN],
            PileupMethod.DEEPTOOLS: [Assay.MCC]
        }

    def _is_compatible(self, method: PileupMethod, scale: DataScalingTechnique, assay: Assay) -> bool:
        if scale in self.incompatible_methods.get(method, []):
            return False
        if assay in self.incompatible_methods.get(method, []):
            return False
        return True

    def generate_bigwig_paths(self) -> list[str]:
        paths = []

        for method in self.pileup_methods:
            for scale in self.scale_methods:
                if not self._is_compatible(method, scale, self.assay):
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
    assay: Assay
    names: list[str]
    peak_calling_method: list[PeakCallingMethod]
    output_dir: str = "seqnado_output"
    is_merged: bool = False

    @property
    def prefix(self) -> str:
        return f"{self.output_dir}/peaks/"

    @field_validator("assay")
    def validate_assay(cls, value):
        if value not in AssaysWithPeakCalling:
            raise ValueError(f"Invalid assay for peak calling: {value}")
        return value

    @property
    def peak_files(self) -> list[str]:
        return expand(
            self.prefix + "{method}/{sample}.bed" if not self.is_merged else self.prefix + "{method}/merged/{sample}.bed",
            sample=self.names,
            method=[m.value for m in self.peak_calling_method],
        )

    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of peak calling files."""
        return self.peak_files


class HeatmapFiles(BaseModel):
    assay: Assay

    @field_validator("assay")
    def validate_assay(cls, value):
        if value not in AssaysWithHeatmaps:
            raise ValueError(f"Invalid assay for heatmap: {value}")
        return value

    output_dir: str = "seqnado_output"

    @property
    def heatmap_files(self) -> list[str]:
        return [
            f"{self.output_dir}/heatmap/heatmap.pdf",
            f"{self.output_dir}/heatmap/metaplot.pdf",
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
    assay: Assay
    names: list[str]
    output_dir: str = "seqnado_output"

    @field_validator("assay")
    def validate_assay(cls, value):
        if value not in AssaysWithSpikein:
            raise ValueError(f"Invalid assay for spike-in: {value}")
        return value

    @property
    def norm_factors(self):
        return f"{self.output_dir}/resources/normalisation_factors.tsv"

    @computed_field
    @property
    def files(self) -> list[str]:
        return [self.norm_factors]


class PlotFiles(BaseModel):
    coordinates: Path
    file_format: Literal["svg", "png", "pdf"] = "svg"
    output_dir: str = "seqnado_output"

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

            outdir = Path(f"{self.output_dir}/genome_browser_plots/")
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
        return [str(p) for p in self.plot_names]

class SNPFilesRaw(BaseModel):
    assay: Assay
    names: list[str]
    output_dir: str = "seqnado_output"

    @property
    def snp_files(self) -> list[str]:
        return expand(
            f"{self.output_dir}/variant/{{sample}}.vcf.gz",
            sample=self.names,
        )
    
    @computed_field
    @property
    def files(self) -> list[str]:
        return self.snp_files

class SNPFilesAnnotated(BaseModel):
    assay: Assay
    names: list[str]
    output_dir: str = "seqnado_output"

    @property
    def anno_snp_files(self) -> list[str]:
        return expand(
            f"{self.output_dir}/variant/{{sample}}.anno.vcf.gz",
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
    output_dir: str = "seqnado_output"

    @property
    def prefix(self) -> str:
        return f"{self.output_dir}/methylation/"

    @property
    def split_bams_files(self) -> List[str]:
        return expand(
            f"{self.output_dir}/aligned/spikein/{{sample}}_{{genome}}.bam",
            sample=self.names,
            genome=self.genomes,
        )
    
    @property
    def methyldackel_files(self) -> List[str]:
        """
        For TAPS, return the inverted file ({sample}_{genome}_CpG_inverted.bedGraph).
        For other methods, return the direct file ({sample}_{genome}_CpG.bedGraph).
        """
        if self.method.value == "taps":
            file_pattern = f"{self.output_dir}/methylation/methyldackel/{{sample}}_{{genome}}_CpG_inverted.bedGraph"
        else:
            file_pattern = f"{self.output_dir}/methylation/methyldackel/{{sample}}_{{genome}}_CpG.bedGraph"
        return expand(
            file_pattern,
            sample=self.names,
            genome=self.genomes,
        )
    
    @property
    def methylation_bias(self) -> List[str]:
        """Return the methylation bias files."""
        files = []
        files.append(f"{self.output_dir}/methylation/methylation_conversion.tsv")
        files.extend(
            expand(
                f"{self.output_dir}/methylation/methyldackel/bias/{{sample}}_{{genome}}.txt",
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
    output_dir: str = "seqnado_output"

    @computed_field
    @property
    def files(self) -> List[str]:
        """Return a list of contact files."""
        # Simplified to just return basic files based on sample names
        # since design_dataframe and viewpoints_grouped are not provided
        return [
            f"{self.output_dir}/contacts/{name}/{name}.mcool" 
            for name in self.names
        ]
    


class QuantificationFiles(BaseModel):
    """Base class for quantification files."""
    
    assay: Assay
    methods: list[QuantificationMethod] = Field(default_factory=list)
    names: list[str]
    groups: SampleGroups
    output_dir: str = "seqnado_output"

    @property
    def prefix(self) -> str:
        return f"{self.output_dir}/quantification"

    @field_validator("methods", mode="before")
    def validate_methods_and_assays(cls, v: list[QuantificationMethod], info) -> list[QuantificationMethod]:
        # In Pydantic v2, use info.data to access other fields
        assay = info.data.get("assay") if hasattr(info, 'data') else None
        if assay == Assay.RNA:
            return [m for m in v if m in [QuantificationMethod.FEATURE_COUNTS, QuantificationMethod.SALMON]]
        return [m for m in v if m == QuantificationMethod.FEATURE_COUNTS]

    @property
    def combined_counts_file(self) -> list[str]:
        """Return the combined read counts file."""
        return expand(self.prefix + "/{methods}/read_counts.tsv", methods=[m.value for m in self.methods])

    @property
    def grouped_counts_files(self) -> list[str]:
        """Return the grouped read counts files."""
        files = []
        for group in self.groups.groups:
            files.extend(
                expand(
                    f"{self.prefix}{group}/read_counts.tsv",
                    group=group,
                    methods=[m.value for m in self.methods],
                )
            )
        return files
    
    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of quantification files."""
        return [*self.combined_counts_file, *self.grouped_counts_files]
            
    

class GeoSubmissionFiles(BaseModel):
    """Class to handle files for GEO submission."""
    
    assay: Assay
    names: list[str]
    seqnado_files: list[str] = Field(default_factory=list, description="Unfiltered list of files. Will be filtered based on allowed extensions and added.")
    allowed_extensions: list[str] = Field(default_factory=lambda: [".bigWig", ".bed", ".tsv", ".vcf.gz"])
    output_dir: str = "seqnado_output"
    
    @property
    def default_files(self):
        return [
            f"{self.output_dir}/geo_submission/md5sums.txt",
            f"{self.output_dir}/geo_submission/raw_data_checksums.txt",
            f"{self.output_dir}/geo_submission/processed_data_checksums.txt",
            f"{self.output_dir}/geo_submission/samples_table.txt",
            f"{self.output_dir}/geo_submission/protocol.txt",
        ]
    
    @property
    def upload_directory(self) -> Path:
        return Path(f"{self.output_dir}/geo_submission") / self.assay.clean_name

    @property
    def upload_instructions(self) -> Path:
        return Path(f"{self.output_dir}/geo_submission") / "upload_instructions.txt"
    

    @property
    def raw_files(self) -> list[str]:
        """Return FASTQ files for raw data."""
        return expand(
            str(self.upload_directory / "{sample}_{read}.fastq.gz"),
            sample=self.names,
            read=["1", "2"],
        )
    
    @property
    def processed_data_files(self) -> list[str]:
        """Return processed files for GEO submission."""
        files = []
        for file in self.seqnado_files:
            file_path = Path(file)
            # Check both single suffix and compound suffix (e.g., .vcf.gz)
            ext = ''.join(file_path.suffixes)
            if ext in self.allowed_extensions or file_path.suffix in self.allowed_extensions:
                # Need to flatten the file un-nest the directory structure
                # e.g. seqnado_output/bigwigs/METHOD/SCALE/NAME.bigWig
                basename = file_path.stem
                scale_method = file_path.parent.name
                method = file_path.parent.parent.name
                
                files.append(
                    str(self.upload_directory / f"{basename}_{method}_{scale_method}{ext}")
                )

        return files
    
    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of all files for GEO submission."""
        return [*self.default_files, *self.raw_files, *self.processed_data_files]
