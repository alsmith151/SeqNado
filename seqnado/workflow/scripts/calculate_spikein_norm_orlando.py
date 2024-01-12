import pandas as pd
import numpy as np
import pathlib
from loguru import logger

# Set up logging
logger.add(snakemake.log[0], level="INFO")

with logger.catch():
    logger.info("Calculating normalization factors")

    # Read in stats
    stats_files = snakemake.input
    all_readcounts = []  

    for stats in stats_files:
        file_path = pathlib.Path(stats)
        readcounts = pd.read_csv(file_path, sep="\t")
        all_readcounts.append(readcounts)

    df_counts = pd.concat(all_readcounts, ignore_index=True)

    # Calculate the ChIP spike-in normalization factor
    df_counts["norm_factor"] = 1 / (df_counts["spikein_reads"] / 1e6)
    # if df_counts["norm_factor"] == inf change to 1
    df_counts["norm_factor"] = df_counts["norm_factor"].replace([np.inf, -np.inf], 1)

    df_counts["scale_factor"] = 1 / df_counts["norm_factor"]
    # if df_counts["scale_factor"] == 0 change to 1
    df_counts["scale_factor"] = df_counts["scale_factor"].replace([0], 1)
    
    df_counts["spikein_percent"] = (df_counts["spikein_reads"] / df_counts["reference_reads"] * 100)
    

    # Save the DataFrame with the calculated normalization factors
    df_counts.to_csv(snakemake.output.normalisation_table, sep="\t", index=False)

    scale = df_counts[["sample", "scale_factor"]].set_index("sample")["scale_factor"]
    scale.to_json(snakemake.output.normalisation_factors)