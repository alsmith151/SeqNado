# import hashlib
# import os
# import pathlib
# import re
# from collections import defaultdict
# from typing import Any, Dict, List, Literal, Optional, Tuple, Union

# import numpy as np
# import pandas as pd
# import snakemake
# from loguru import logger
# from pydantic import BaseModel, Field, computed_field, validator
# from pydantic.dataclasses import dataclass
# from snakemake.io import expand







# def define_output_files(
#     snakemake_design: Union[Design, DesignIP],
#     assay: Literal["ChIP", "ATAC", "RNA", "SNP"],
#     fastq_screen: bool = False,
#     chip_spikein_normalisation: bool = False,
#     sample_names: list = None,
#     pileup_method: list = None,
#     peak_calling_method: list = None,
#     make_bigwigs: bool = False,
#     call_peaks: bool = False,
#     make_heatmaps: bool = False,
#     make_ucsc_hub: bool = False,
#     call_snps: bool = False,
#     annotate_snps: bool = False,
#     **kwargs,
# ) -> list:
#     """Define output files for the pipeline"""

#     analysis_output = [
#         "seqnado_output/qc/fastq_raw_qc.html",
#         "seqnado_output/qc/fastq_trimmed_qc.html",
#         "seqnado_output/qc/alignment_raw_qc.html",
#         "seqnado_output/qc/alignment_filtered_qc.html",
#         "seqnado_output/design.csv",
#     ]
#     assay_output = []

#     if fastq_screen:
#         assay_output.append("seqnado_output/qc/full_fastqscreen_report.html")

#     if kwargs["remove_pcr_duplicates_method"] == "picard":
#         analysis_output.append("seqnado_output/qc/library_complexity_qc.html")

#     if make_heatmaps:
#         assay_output.extend(
#             [
#                 "seqnado_output/heatmap/heatmap.pdf",
#                 "seqnado_output/heatmap/metaplot.pdf",
#             ]
#         )

#     if make_ucsc_hub:
#         hub_dir = pathlib.Path(kwargs["ucsc_hub_details"]["directory"])
#         hub_name = kwargs["ucsc_hub_details"]["name"]
#         hub_txt = hub_dir / f"{hub_name}.hub.txt"
#         analysis_output.append(str(hub_txt))

#     if assay in ["ChIP", "ATAC"]:
#         if assay == "ChIP" and chip_spikein_normalisation:
#             if assay == "ChIP":
#                 assay_output.extend(
#                     [
#                         "seqnado_output/qc/full_fastqscreen_report.html",
#                         "seqnado_output/resources/normalisation_factors.tsv",
#                     ]
#                 )

#         if make_bigwigs and pileup_method:
#             assay_output.extend(
#                 expand(
#                     "seqnado_output/bigwigs/{method}/{sample}.bigWig",
#                     sample=sample_names,
#                     method=pileup_method,
#                 )
#             )

#         if call_peaks and peak_calling_method:
#             if assay == "ChIP":
#                 # Add peak calling output
#                 assay_output.extend(
#                     expand(
#                         "seqnado_output/peaks/{method}/{ip}.bed",
#                         ip=kwargs["ip"],
#                         method=peak_calling_method,
#                     )
#                 )

#             else:
#                 assay_output.extend(
#                     expand(
#                         "seqnado_output/peaks/{method}/{sample}.bed",
#                         sample=sample_names,
#                         method=peak_calling_method,
#                     )
#                 )

#         if "merge" in snakemake_design.to_dataframe().columns:
#             for group_name, df in snakemake_design.to_dataframe().groupby("merge"):
#                 assay_output.extend(
#                     expand(
#                         "seqnado_output/peaks/consensus/{group_name}.bed",
#                         group_name=group_name,
#                     )
#                 )

#     elif assay == "RNA":
#         assay_output.extend(
#             [
#                 "seqnado_output/feature_counts/read_counts.tsv",
#                 *expand(
#                     "seqnado_output/aligned/{sample}.bam",
#                     sample=sample_names,
#                 ),
#             ]
#         )

#         if make_bigwigs and pileup_method:
#             assay_output.extend(
#                 expand(
#                     "seqnado_output/bigwigs/{method}/{sample}_{strand}.bigWig",
#                     sample=sample_names,
#                     method=pileup_method,
#                     strand=["plus", "minus"],
#                 )
#             )

#         if kwargs["run_deseq2"]:
#             if kwargs["can_run_deseq2"]:
#                 assay_output.append(f"deseq2_{kwargs['project_name']}.html")
#             else:
#                 logger.warning(
#                     "Not running DESeq2 as no 'deseq2' column in design file."
#                 )

#     elif assay == "SNP":
#         if call_snps:
#             assay_output.expand(
#                 "seqnado_output/variant/{sample}_filtered.anno.vcf.gz",
#                 sample=sample_names,
#             )

#         if annotate_snps:
#             assay_output.append(
#                 expand(
#                     "seqnado_output/variant/{sample}_filtered.stats.txt",
#                     sample=sample_names,
#                 ),
#             )

#     analysis_output.extend(assay_output)

#     return analysis_output
