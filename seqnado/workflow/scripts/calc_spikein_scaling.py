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



# List all .txt files from split_dir
files = [f for f in os.listdir(split_dir) if f.endswith(".txt")]
files.sort()

dfs = []
for report in files:
    file_path = os.path.join(split_dir, report)
    df = pd.read_csv(file_path, sep="\t")
    if (
        "Sample" not in df.columns
        or "n_sample" not in df.columns
        or "n_exogenous" not in df.columns
    ):
        ConnectionRefusedError
    dfs.append(df)

if dfs:
    scale_factors = pd.concat(dfs, ignore_index=True)
    scale_factors.to_csv(
        os.path.join(split_dir, "scale_factors.csv"),
        index=False,
    )
else:
    print("No valid data found in the input files.")

# subset columns of interest
scale_factors = scale_factors[["Sample", "n_sample", "n_exogenous"]]

# Calculate the ChIP spike-in normalization factor
scale_factors["chip_spike_in_norm_factor"] = 1 / (scale_factors["n_exogenous"] / 1e6)
scale_factors["exo_perc"] = (
    scale_factors["n_exogenous"] / scale_factors["n_sample"] * 100
)

# Save the DataFrame with the calculated normalization factors
scale_factors.to_csv(
    os.path.join(split_dir, "scale_factors_with_norm_factors.csv"),
    index=False,
)

scale_factors.to_csv(
    os.path.join(main_dir, "scale_factors.csv"),
    index=False,
)
