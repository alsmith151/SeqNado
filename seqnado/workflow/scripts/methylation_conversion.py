import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from loguru import logger

# Handle Snakemake execution context
try:
    input_files = snakemake.input
    output_file = snakemake.output[0]
    assay = snakemake.params.assay
    log_file = snakemake.log[0]
    logger.add(log_file, level="INFO")
except NameError:
    input_files = sys.argv[1:-3]  # Capture multiple input files
    output_file = sys.argv[-3]
    assay = sys.argv[-2]
    log_file = sys.argv[-1]
    logger.add(log_file, level="INFO")


def plot_methylation_bias(summary_df, output_file):
    """
    Generates a bar plot for methylation bias summary with a merged legend placed below the subplots.

    Parameters:
        summary_df (pd.DataFrame): DataFrame containing conversion rates.
        output_file (str): Path to save the output plot.
    """
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
    print(f"Methylation bias summary plot saved: {output_file}")


def combine_methylation_bias(input_files, output_file, assay):
    """
    Combines methylation bias statistics from multiple samples into a single summary file.

    Parameters:
        input_files (list): List of input file paths containing bias statistics.
        output_file (str): Path to the output summary file.
        assay (str): Type of assay ("bisulfite" or "taps").
    """
    combined_data = []

    for bias_file in input_files:
        df = pd.read_csv(bias_file, sep="\t")

        # Extract sample & genome from filename
        filename = os.path.basename(bias_file).replace(".txt", "")
        parts = filename.split("_")
        sample = parts[0]
        genome = "_".join(parts[1:])  # Handle cases where genome has "_"

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

        # Adjust for TAPS assay
        if assay == "taps":
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


if __name__ == "__main__":
    combine_methylation_bias(input_files, output_file, assay)
