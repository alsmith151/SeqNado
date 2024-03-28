import pandas as pd
import pathlib
from loguru import logger

# Set up logging
logger.add(snakemake.log[0], level="INFO")

with logger.catch():
    logger.info("Collating read counts from Salmon quantification")

    files = snakemake.input.counts
    all_readcounts = []

    for file in files:
        file_path = pathlib.Path(file)
        sample_name = file_path.parts[-2].replace("salmon_", "")
        readcounts = pd.read_csv(file_path, sep="\t")
        readcounts.rename(columns={"NumReads": sample_name}, inplace=True)
        readcounts = readcounts[["Name", sample_name]]
        all_readcounts.append(readcounts)

    cumulative_df = all_readcounts[0]
    for readcounts in all_readcounts[1:]:
        cumulative_df = pd.merge(
            cumulative_df,
            readcounts,
            on=["Name"],
            how="outer",
        )

    cumulative_df.to_csv(snakemake.output.count_table, sep="\t", index=False)
