import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from loguru import logger


def snakemake_setup():
    if "snakemake" not in globals():
        raise RuntimeError("This script must be run via Snakemake.")
    log_file = snakemake.log[0]
    logger.remove()
    logger.add(log_file, format="{time} {level} {message}", level="INFO")
    logger.add(sys.stderr, format="{time} {level} {message}", level="ERROR")
    input_files = snakemake.input
    output_file = snakemake.output[0]
    method = snakemake.params.method
    return input_files, output_file, method


def plot_methylation_bias(summary_df, output_file):
    """
    Generates a bar plot for methylation bias summary with a merged legend placed below the subplots.

    Parameters:
        summary_df (pd.DataFrame): DataFrame containing conversion rates.
        output_file (str): Path to save the output plot.
    """
    logger.info(f"Generating methylation bias plot: {output_file}")
    sns.set_theme(style="whitegrid", palette="Set2")
    fig, ax = plt.subplots(1, 2, figsize=(12, 7), sharey=True)
    sns.barplot(
        data=summary_df,
        x="Sample",
        y="Read1_Conversion",
        hue="Genome",
        ax=ax[0],
    )
    ax[0].set_title("Read 1", fontsize=14, fontweight="bold")
    ax[0].set_ylabel("Conversion Rate (%)", fontsize=12)
    ax[0].set_xlabel("")
    ax[0].set_ylim(0, 100)

    sns.barplot(
        data=summary_df,
        x="Sample",
        y="Read2_Conversion",
        hue="Genome",
        ax=ax[1],
    )
    ax[1].set_title("Read 2", fontsize=14, fontweight="bold")
    ax[1].set_ylabel("Conversion Rate (%)", fontsize=12)
    ax[1].set_xlabel("")
    ax[1].set_ylim(0, 100)

    for a in ax:
        a.set_xticklabels(a.get_xticklabels(), rotation=45, ha="right")

    handles, labels = ax[0].get_legend_handles_labels()
    legend = fig.legend(
        handles,
        labels,
        title="Genome",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=len(labels),
        fontsize=12,
        title_fontsize=12,
        frameon=False,
    )

    ax[0].legend_.remove()
    ax[1].legend_.remove()

    plt.tight_layout()
    plt.savefig(output_file, bbox_extra_artists=(legend,), bbox_inches="tight")
    plt.close()
    logger.info(f"Methylation bias summary plot saved: {output_file}")


def combine_methylation_bias(input_files: list, output_file: str, method: str):
    """
    Combines methylation bias statistics from multiple samples into a single summary file.

    Parameters:
        input_files (list): List of input file paths containing bias statistics.
        output_file (str): Path to the output summary file.
        method (METHOD): Type of method ("bisulfite" or "taps").
    """
    logger.info(f"Combining methylation bias from {len(input_files)} input files.")
    combined_data = []

    for bias_file in input_files:
        logger.info(f"Processing bias file: {bias_file}")
        try:
            df = pd.read_csv(bias_file, sep="\t")
        except Exception as e:
            logger.error(f"Failed to read {bias_file}: {e}")
            continue

        # Extract sample & genome from filename
        filename = os.path.basename(bias_file).replace(".txt", "")
        parts = filename.split("_")
        sample = parts[0]
        genome = "_".join(parts[1:])  # Handle cases where genome has "_"
        logger.debug(f"Sample: {sample}, Genome: {genome}")

        # Compute aggregated methylation bias for Read 1 and Read 2
        read1 = df[df["Read"] == 1].sum()
        read2 = df[df["Read"] == 2].sum()

        # Calculate conversion rate
        read1_meth = (
            (read1["nMethylated"] / (read1["nMethylated"] + read1["nUnmethylated"]))
            * 100
            if (read1["nMethylated"] + read1["nUnmethylated"]) > 0
            else None
        )
        read2_meth = (
            (read2["nMethylated"] / (read2["nMethylated"] + read2["nUnmethylated"]))
            * 100
            if (read2["nMethylated"] + read2["nUnmethylated"]) > 0
            else None
        )

        # Adjust for TAPS method
        if method == "taps":
            read1_meth = 100 - read1_meth if read1_meth is not None else None
            read2_meth = 100 - read2_meth if read2_meth is not None else None

        # Store result
        combined_data.append(
            [
                sample,
                genome,
                read1["nMethylated"],
                read1["nUnmethylated"],
                read1_meth,
                read2["nMethylated"],
                read2["nUnmethylated"],
                read2_meth,
            ]
        )

    # Convert to DataFrame
    summary_df = pd.DataFrame(
        combined_data,
        columns=[
            "Sample",
            "Genome",
            "Read1_nMethylated",
            "Read1_nUnmethylated",
            "Read1_Conversion",
            "Read2_nMethylated",
            "Read2_nUnmethylated",
            "Read2_Conversion",
        ],
    )

    # Save as tab-separated file
    summary_df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"Methylation bias summary saved: {output_file}")

    # Convert "None" to NaN before filtering
    summary_df.replace({None: pd.NA}, inplace=True)
    summary_df.dropna(subset=["Read1_Conversion", "Read2_Conversion"], inplace=True)
    try:
        plot_file = output_file.replace(".tsv", ".png")
        plot_methylation_bias(summary_df, plot_file)
    except Exception as e:
        logger.error(f"Error generating plot: {e}")


def main():
    input_files, output_file, method = snakemake_setup()
    logger.info(
        f"Starting methylation conversion script with input_files={input_files}, output_file={output_file}, method={method}"
    )
    combine_methylation_bias(input_files, output_file, method)

if __name__ == "__main__":
    main()
