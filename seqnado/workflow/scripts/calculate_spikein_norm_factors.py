import os
import sys
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

# Calculate readcounts for reference and spikein samples
bam_ref = [pathlib.Path(p) for p in snakemake.input.bam_ref]
bam_spikein = [pathlib.Path(p) for p in snakemake.input.bam_spikein]
readcounts_ref = get_readcounts(bam_ref)
readcounts_spikein = get_readcounts(bam_spikein)

# Read in design matrix
df_design_raw = (pd.read_csv(snakemake.input.design, sep=",")
               .assign(ip=lambda df: df["sample"] + "_" + df["antibody"])
               [["sample", "antibody", "ip", "control"]]
            )

df_design = df_design_raw.melt(id_vars=["sample", "antibody"], var_name="type", value_name="name")

# Merge readcounts with design matrix
df_counts = (pd.DataFrame([readcounts_ref, readcounts_spikein], index=["ref", "spikein"])
               .T
               )

df_counts = (df_design.merge(df_counts, left_on="name", right_index=True))

# Pivot for easier handling
df_counts = df_counts.pivot_table(index=["sample", "antibody"], columns="type", values=["ref", "spikein"])
df_counts.columns = ["_".join(col) for col in df_counts.columns.values]
df_counts = df_counts.reset_index()


# Calculate normalization factors
df_counts = df_counts.assign(
    reads_per_spikein_ip=lambda df: df["ref_ip"] / df["spikein_ip"],
    reads_per_spikein_control=lambda df: df["ref_control"] / df["spikein_control"],
    relative_signal=lambda df: df["reads_per_spikein_ip"] / df["reads_per_spikein_control"],
    reads_per_million_ip=lambda df: 10e6 / df["ref_ip"],
    reads_per_million_control=lambda df: 10e6 / df["ref_control"],
    norm_factor=lambda df: df["relative_signal"] * df["reads_per_million_ip"]
)

# Write out normalization factors
df_counts = df_counts.merge(df_design_raw, on=["sample", "antibody"])
df_counts.to_csv(snakemake.output.normalisation_table, sep="\t", index=False)

norm_ip = df_counts[["ip", "norm_factor"]].set_index("ip")["norm_factor"]
norm_control = df_counts[["control", "norm_factor"]].set_index("control").rename(columns={"reads_per_million_control": "norm_factor"})["norm_factor"]

norm = pd.concat([norm_ip, norm_control])
norm.to_json(snakemake.output.normalisation_factors)








