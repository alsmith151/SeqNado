import pandas as pd
import pathlib
import pysam
from typing import List
from loguru import logger

# Set up logging
logger.add(snakemake.log[0], level="INFO")


def get_readcounts(bam_files: List[pathlib.Path]):
    readcounts = {}
    for bam_file in bam_files:
        bam = pysam.AlignmentFile(bam_file, "rb")
        readcounts[bam_file.stem] = bam.mapped
    return pd.Series(readcounts)


with logger.catch():
    logger.info("Calculating normalization factors")

    # Calculate readcounts for reference and spikein samples
    bam_ref = [pathlib.Path(p) for p in snakemake.input.bam_ref]
    bam_spikein = [pathlib.Path(p) for p in snakemake.input.bam_spikein]
    readcounts_ref = get_readcounts(bam_ref)
    readcounts_spikein = get_readcounts(bam_spikein)

    # Read in design matrix
    df_design = pd.read_csv(snakemake.input.design, sep=",")[["sample", "antibody"]]

    # Merge readcounts with design matrix
    df_counts = pd.DataFrame(
        [readcounts_ref, readcounts_spikein], index=["reference", "spikein"]
    )

    df_counts = df_design.merge(df_counts, left_on="name", right_index=True)

    # Pivot for easier handling
    df_counts = df_counts.pivot_table(
        index=["sample"], columns="type", values=["reference", "spikein"]
    )
    df_counts.columns = ["_".join(col) for col in df_counts.columns.values]
    df_counts = df_counts.reset_index()

    # subset columns of interest
    df_counts = df_counts[["sample", "reference", "spikein"]]

    df_counts["chip_spike_in_norm_factor"] = 1 / (df_counts["spikein"] / 1e6)
    df_counts["exo_perc"] = (
        df_counts["spikein"] / df_counts["reference"] * 100
    )

    # Write out normalization factors
    df_counts.to_csv(snakemake.output.normalisation_table, sep="\t", index=False)
