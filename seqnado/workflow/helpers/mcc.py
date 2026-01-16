"""Helper functions for MCC (Multi-way Chromatin Contact) workflows."""

import json
from pathlib import Path
from typing import Dict, List
from itertools import chain, product, combinations
from snakemake.io import expand


def get_n_cis_scaling_factor(wc, OUTPUT_DIR):
    """
    Calculate the n_cis scaling factor for MCC normalization.

    Args:
        wc: Snakemake wildcards object containing 'viewpoint_group' and optionally 'group' or 'sample'.
        OUTPUT_DIR: The output directory path.

    Returns:
        float: The n_cis scaling factor, or 1 if stats file not found, or 0 if n_total is 0.
    """
    # Can either extract the stats for the sample or the group
    if hasattr(wc, "group"):
        stats_file = OUTPUT_DIR + f"/resources/{wc.group}_ligation_stats.json"
    else:
        # If not, then use the sample
        stats_file = OUTPUT_DIR + f"/resources/{wc.sample}_ligation_stats.json"

    # Create Path object and ensure the file exists
    stats_path = Path(stats_file)
    if not stats_path.exists():
        return 1

    with open(stats_path, "r") as r:
        stats = json.load(r)

    # Check if viewpoint_group exists in stats
    if wc.viewpoint_group not in stats:
        raise KeyError(
            f"Viewpoint group '{wc.viewpoint_group}' not found in stats file"
        )

    # Check if required keys exist in the viewpoint group stats
    required_keys = ["n_cis", "n_total"]
    missing_keys = [
        key for key in required_keys if key not in stats[wc.viewpoint_group]
    ]
    if missing_keys:
        raise KeyError(f"Missing required keys in stats: {', '.join(missing_keys)}")

    # Avoid division by zero
    if stats[wc.viewpoint_group]["n_total"] == 0:
        return 0

    return (
        stats[wc.viewpoint_group]["n_cis"] / stats[wc.viewpoint_group]["n_total"]
    ) * 1e6


def get_mcc_bam_files_for_merge(wildcards, SAMPLE_GROUPINGS, OUTPUT_DIR):
    """
    Get BAM files for merging based on sample names.

    Args:
        wildcards: Snakemake wildcards object containing 'group'.
        SAMPLE_GROUPINGS: The sample groupings object.
        OUTPUT_DIR: The output directory path.

    Returns:
        list: List of BAM file paths to merge, or empty list if group not found.
    """
    try:
        group = SAMPLE_GROUPINGS.get_grouping("consensus").get_group(wildcards.group)
        sample_names = group.samples
        bam_files = [
            OUTPUT_DIR + f"/mcc/replicates/{sample}/{sample}.bam"
            for sample in sample_names
        ]
        return bam_files
    except KeyError:
        return []


def identify_extracted_bam_files(wildcards, checkpoints):
    """
    Identify extracted BAM files from checkpoint output.

    Args:
        wildcards: Snakemake wildcards object.
        checkpoints: Snakemake checkpoints object.

    Returns:
        list: List of extracted BAM file paths.
    """
    checkpoint_output = checkpoints.identify_viewpoint_reads.get(**wildcards)
    outdir = Path(checkpoint_output.output.bams)

    from snakemake.io import glob_wildcards

    viewpoints = glob_wildcards(str(outdir / "{viewpoint}.bam")).viewpoint

    return [str(outdir / f"{viewpoint}.bam") for viewpoint in viewpoints]


def redefine_viewpoints(samples, checkpoints):
    """
    Redefine the set of viewpoints to be the intersection of viewpoints across all samples.

    The issue is that some viewpoints may not be present in all samples or may not have enough reads to be considered.

    Args:
        samples: List of sample names.
        checkpoints: Snakemake checkpoints object.

    Returns:
        list: List of viewpoints common to all samples.
    """
    from snakemake.io import glob_wildcards

    viewpoint_set = set()

    for ii, sample in enumerate(samples):
        checkpoint_output = checkpoints.identify_viewpoint_reads.get(sample=sample)
        outdir = Path(checkpoint_output.output.bams)
        viewpoints = glob_wildcards(str(outdir / "{viewpoint}.bam")).viewpoint

        if ii == 0:
            viewpoint_set = set(viewpoints)
        else:
            viewpoint_set = viewpoint_set.intersection(viewpoints)

    return list(viewpoint_set)


def extract_viewpoints(viewpoints_path: str) -> List[str]:
    """
    Extracts the viewpoints from the config.
    """
    import numpy as np
    import pandas as pd

    # Read BED file using pandas (pyranges replacement)
    bed_columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
    try:
        df = pd.read_csv(viewpoints_path, sep="\t", header=None, comment="#")
        # Assign column names based on the number of columns
        df.columns = bed_columns[: len(df.columns)]
        # Ensure we have at least the minimum required columns
        if "Name" not in df.columns:
            df["Name"] = (
                df["Chromosome"]
                + ":"
                + df["Start"].astype(str)
                + "-"
                + df["End"].astype(str)
            )
    except Exception as e:
        raise ValueError(f"Error reading BED file {viewpoints_path}: {e}")

    df = df.assign(
        viewpoint=lambda df: np.where(
            df["Name"].str.contains(r"-chr.*?-\d+-\d+$"),
            df["Name"],
            df["Name"]
            + "-"
            + df["Chromosome"].astype(str)
            + "-"
            + df["Start"].astype(str)
            + "-"
            + df["End"].astype(str),
        )
    )

    viewpoints = set(df["viewpoint"].tolist())
    return viewpoints


def viewpoint_to_grouped_viewpoint(viewpoints: List[str]) -> Dict[str, str]:
    """
    Groups the viewpoints that consist of multiple oligos into a dictionary.
    """
    import re

    has_coordinate = re.compile(r"(.*?)-chr([0-9]+|X|Y|M|MT)-\d+-\d+$")
    viewpoint_to_grouped_mapping = {}

    for viewpoint in viewpoints:
        has_coordinate_match = has_coordinate.match(viewpoint)
        if has_coordinate_match:
            viewpoint_name = has_coordinate_match.group(1)
            viewpoint_to_grouped_mapping[viewpoint] = viewpoint_name

    return viewpoint_to_grouped_mapping

