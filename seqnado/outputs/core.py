from itertools import chain
from pathlib import Path
from typing import Any, List, Literal, Union

import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field, computed_field
from snakemake.io import expand

from seqnado import (
    Assay,
)
from seqnado.config import SeqnadoConfig
from seqnado.inputs import IPSampleCollection, SampleCollection, SampleGroups
from seqnado.outputs.files import (
    BigWigFiles,
    FileCollection,
    HeatmapFiles,
    HubFiles,
    PeakCallingFiles,
    PlotFiles,
    QCFiles,
    QuantificationFiles,
    SpikeInFiles,
    SNPFilesRaw,
    SNPFilesAnnotated,
    MethylationFiles,
    ContactFiles,
)


class OutputFiles(BaseModel):
    """A collection of output files for a SeqNado project."""

    files: list[str] = Field(default_factory=list)

    @property
    def all_files(self) -> List[str]:
        """Return all files in the output collection."""
        return self.files

    @property
    def bigwig_files(self) -> List[str]:
        """Return only bigwig files."""
        return [f for f in self.files if f.endswith(".bigWig")]

    @property
    def peak_files(self) -> List[str]:
        """Return only peak calling files."""
        return [f for f in self.files if f.endswith(".bed")]

    @property
    def heatmap_files(self) -> List[str]:
        """Return only heatmap files."""
        return [f for f in self.files if f.endswith(".pdf") and "heatmap" in f]


class OutputBuilder:
    def __init__(
        self,
        assay: Assay,
        samples: SampleCollection | IPSampleCollection,
        config: SeqnadoConfig,
        sample_groups: SampleGroups | None = None,
    ):
        """Initializes the OutputBuilder with the given assay, samples, and configuration.
        Args:
            assay (Assay): The type of assay being processed.
            samples (SampleCollection | IPSampleCollection): The collection of samples to process.
            config (SeqnadoConfig): The configuration for the SeqNado project.
            sample_groups (Optional[SampleGroups]): Optional groups of samples for merging.
                If provided, will be used to create grouped bigwig files and grouped peak files.
        Raises:
            ValueError: If the provided assay is not supported or if sample groups are provided
                but not defined in the configuration.
        """

        self.assay = assay
        self.samples = samples
        self.config = config
        self.sample_groups = sample_groups

        # Initialize an empty list to hold file collections
        self.file_collections: list[FileCollection] = []

    def add_qc_files(self) -> None:
        qc_files = QCFiles(assay=self.assay, samples=self.samples)
        self.file_collections.append(qc_files)

    def add_individual_bigwig_files(self) -> None:
        bigwig_files = BigWigFiles(
            assay=self.assay,
            names=self.samples.sample_names,
            pileup_methods=self.config.assay_config.bigwigs.pileup_method,
            scale_methods=self.config.assay_config.bigwigs.scale_methods,
        )
        self.file_collections.append(bigwig_files)

    def add_grouped_bigwig_files(self) -> None:
        for group in self.sample_groups.groups:
            bigwig_files = BigWigFiles(
                assay=self.assay,
                names=[group.name],
                pileup_methods=self.config.assay_config.bigwigs.pileup_method,
                scale_methods=self.config.assay_config.bigwigs.scale_methods,
            )
            self.file_collections.append(bigwig_files)

    def add_peak_files(self) -> None:
        """Add peak files to the output collection."""
        peaks = PeakCallingFiles(
            assay=self.assay,
            names=self.samples.sample_names,
            peak_calling_method=self.config.assay_config.peak_calling_methods,
        )
        self.file_collections.append(peaks)
    
    def add_grouped_peak_files(self) -> None:
        """Add grouped peak files to the output collection."""
        if not self.sample_groups:
            raise ValueError("Sample groups must be defined to create grouped peak files.")
        
        for group in self.sample_groups.groups:
            peaks = PeakCallingFiles(
                assay=self.assay,
                names=[group.name],
                peak_calling_method=self.config.assay_config.peak_calling_methods,
            )
            self.file_collections.append(peaks)
    
    def add_bigbed_files(self) -> None:
        """Add bigBed files to the output collection."""
        bigbed_files = BigWigFiles(
            assay=self.assay,
            names=self.samples.sample_names,
            pileup_methods=self.config.assay_config.bigwigs.pileup_method,
            scale_methods=self.config.assay_config.bigwigs.scale_methods,
            prefix=Path("seqnado_output/bigbeds/"),
        )
        self.file_collections.append(bigbed_files)

    def add_heatmap_files(self) -> None:
        """Add heatmap files to the output collection."""
        heatmaps = HeatmapFiles(assay=self.assay)
        self.file_collections.append(heatmaps)

    def add_hub_files(self) -> None:
        """Add hub files to the output collection."""
        hub_files = HubFiles(
            hub_dir=Path("seqnado_output/hub"),
            hub_name=self.config.assay_config.ucsc_hub.name,
        )
        self.file_collections.append(hub_files)

    def add_spikein_files(self) -> None:
        """Add spike-in files to the output collection."""
        spikein_files = SpikeInFiles(
            assay=self.assay,
            names=self.samples.sample_names,
        )
        self.file_collections.append(spikein_files)

    def add_plot_files(
        self, coordinates: Path, file_format: Literal["svg", "png", "pdf"] = "svg"
    ) -> None:
        """Add plot files to the output collection."""
        plot_files = PlotFiles(
            coordinates=coordinates,
            file_format=file_format,
        )
        self.file_collections.append(plot_files)
    
    def add_snp_files(self) -> None:
        """Add SNP files to the output collection."""
        snp_files_raw = SNPFilesRaw(
            assay=self.assay,
            names=self.samples.sample_names,
        )
        self.file_collections.append(snp_files_raw)

        if self.config.assay_config.snp_calling.annotate_snps:
            snp_files_annotated = SNPFilesAnnotated(
                assay=self.assay,
                names=self.samples.sample_names,
            )
            self.file_collections.append(snp_files_annotated)
    
    def add_methylation_files(self) -> None:
        """Add methylation files to the output collection."""
        methylation_files = MethylationFiles(
            assay=self.assay,
            names=self.samples.sample_names,
            method=self.config.assay_config.methylation.method,
        )
        self.file_collections.append(methylation_files)
    
    def add_contact_files(self) -> None:
        """Add contact files to the output collection."""
        contact_files = ContactFiles(
            assay=self.assay,
            names=self.samples.sample_names,
            prefix=Path("seqnado_output/mcc/"),
        )
        self.file_collections.append(contact_files)
    
    def add_quantification_files(self) -> None:
        """Add quantification files to the output collection."""
        quantification_files = QuantificationFiles(
            assay=self.assay,
            methods=self.config.assay_config.quantification.methods,
            names=self.samples.sample_names,
            groups=self.sample_groups,
            prefix=Path("seqnado_output/quantification/"),
        )
        self.file_collections.append(quantification_files)
    

    def build(self) -> OutputFiles:
        """Builds the output files collection based on the added file collections."""
        all_files = list(
            chain.from_iterable(p.files for p in self.file_collections if p.files)
        )
        return OutputFiles(files=all_files)










# def to_geo_dataframe(
#     self,
#     assay: Literal["ATAC", "RNA", "SNP"],
#     pipeline_config: dict[str, Any],
# ) -> pd.DataFrame:
#     """
#     Generate GEO-compliant metadata table for submission.

#     Args:
#         assay: One of "ATAC", "RNA", "SNP".
#         pipeline_config: Must contain genome.name and optionally instrument_model.

#     Returns:
#         pandas DataFrame ready for GEO.
#     """
#     # Build GEOSample entries
#     geo_samples: list[GEOSample] = []
#     df = self.to_dataframe()
#     for row in df.itertuples(index=False):
#         # Determine processed files
#         if assay == "RNA":
#             proc_files = [
#                 "read_counts.tsv",
#                 f"{row.sample_name}_plus.bw",
#                 f"{row.sample_name}_minus.bw",
#             ]
#         else:
#             proc_files = [f"{row.sample_name}.bw"]

#         raw_files = [pathlib.Path(row.r1).name]
#         if row.r2:
#             raw_files.append(pathlib.Path(row.r2).name)

#         sample = GEOSample(
#             assay=assay,
#             library_name=row.sample_name,
#             title=row.sample_name,
#             organism=predict_organism(pipeline_config["genome"]["name"]),
#             cell_line=None,
#             cell_type=None,
#             antibody=None,
#             genotype=None,
#             treatment=getattr(row, "treatment", None),
#             time=getattr(row, "time", None),
#             single_or_paired=("paired-end" if row.r2 else "single"),
#             instrument_model=pipeline_config.get(
#                 "instrument_model", "Illumina NovaSeq X"
#             ),
#             description=None,
#             processed_data_file=proc_files,
#             raw_file=raw_files,
#         )
#         geo_samples.append(sample)

#     return GEOSamples(samples=geo_samples).to_dataframe()

# def to_geo_dataframe(
#     self,
#     assay: Literal["ChIP", "CAT"],
#     pipeline_config: dict[str, Any],
# ) -> pd.DataFrame:
#     """
#     Generate GEO metadata for IP/control experiments.

#     Each sample (IP and control) becomes a row.
#     """
#     geo_samples: list[GEOSample] = []
#     df = self.to_dataframe()
#     for rec in df.itertuples(index=False):
#         for side in ("ip", "control"):  # produce both
#             label = getattr(rec, side)
#             if label is None:
#                 continue

#             r1 = getattr(rec, f"{side}_r1")
#             r2 = getattr(rec, f"{side}_r2")
#             raw_files = [pathlib.Path(r1).name]
#             if r2:
#                 raw_files.append(pathlib.Path(r2).name)

#             proc_file = f"{rec.sample_name}_{label}.bigWig"
#             sample = GEOSample(
#                 assay=assay,
#                 library_name=rec.sample_name,
#                 title=rec.sample_name,
#                 organism=predict_organism(pipeline_config["genome"]["name"]),
#                 cell_line=None,
#                 cell_type=None,
#                 antibody=label,
#                 genotype=None,
#                 treatment=getattr(rec, "treatment", None),
#                 time=getattr(rec, "time", None),
#                 single_or_paired=("paired-end" if r2 else "single"),
#                 instrument_model=pipeline_config.get(
#                     "instrument_model", "Illumina NovaSeq X"
#                 ),
#                 description=None,
#                 processed_data_file=[proc_file],
#                 raw_file=raw_files,
#             )
#             geo_samples.append(sample)

#     return GEOSamples(samples=geo_samples).to_dataframe()
