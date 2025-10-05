import pandas as pd
from pathlib import Path
from loguru import logger

# Set up logging
logger.add(snakemake.log[0], level="INFO")

with logger.catch():
    logger.info("Calculating normalization factors")

# Read in stats
stats_files = snakemake.input

all_readcounts = []

for stats in stats_files:
    file_path = Path(stats)
    readcounts = pd.read_csv(file_path, sep="\t")
    all_readcounts.append(readcounts)


df_counts = pd.concat(all_readcounts, ignore_index=True)
df_counts["sample_name"] = df_counts["sample"].str.split("_", expand=True)[0]
df_counts["ip"] = df_counts["sample"].str.split("_", expand=True)[1]
df_counts.drop(columns=["sample"], inplace=True)
df_counts['ip'] = df_counts['ip'].apply(lambda x: x.lower() if x == 'Input' else x)
df_counts = df_counts[['sample_name', 'ip', 'reference_reads', 'spikein_reads']]

df_counts_ip = df_counts[df_counts['ip'] != 'input']
df_counts_input = df_counts[df_counts['ip'] == 'input'].drop(columns=['ip'])
df_counts = df_counts_ip.merge(df_counts_input, on="sample_name", suffixes=('_ip', '_control'))

# Calculate normalization factors
df_counts = df_counts.assign(
    reads_per_spikein_ip=lambda df: df["reference_reads_ip"] / df["spikein_reads_ip"],
    reads_per_spikein_control=lambda df: df["reference_reads_control"] / df["spikein_reads_control"],
    relative_signal=lambda df: df["reads_per_spikein_ip"]
    / df["reads_per_spikein_control"],
    reads_per_million_ip=lambda df: 10e6 / df["reference_reads_ip"],
    reads_per_million_control=lambda df: 10e6 / df["reference_reads_control"],
    norm_factor=lambda df: df["relative_signal"] * df["reads_per_million_ip"],
)


df_counts['sample'] = df_counts['sample_name'] + "_" + df_counts['ip']
# Write out normalization factors
df_counts.to_csv(snakemake.output.normalisation_table, sep="\t", index=False)

norm_ip = df_counts[["sample", "norm_factor"]].set_index(["sample"])["norm_factor"]
df_counts_input['sample'] = df_counts_input['sample_name'] + "_Input"
df_counts_input['norm_factor'] = 1
norm_input = df_counts_input[["sample", "norm_factor"]].set_index(["sample"])["norm_factor"]
norm = pd.concat([norm_ip, norm_input], axis=0)
# unique norm 
norm = norm[~norm.index.duplicated(keep='first')].sort_index(
    ascending=True
)
norm.to_json(snakemake.output.normalisation_factors)