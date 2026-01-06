"""
Create a design matrix for MAGeCK MLE analysis from sample metadata.

This script generates a design matrix by one-hot encoding categorical variables
from the sample metadata. The design matrix is used by MAGeCK MLE to model
complex experimental designs with multiple conditions.

For CRISPR screens, common experimental factors include:
- treatment (e.g., treated vs control)
- timepoint (e.g., day0, day7, day14)
- condition (e.g., drug vs vehicle)
- replicate groups

The script uses the 'deseq2' metadata column if available, which commonly
contains experimental group annotations.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from loguru import logger
import sys

# Set up logging
logger.remove()  # Remove default handler
logger.add(snakemake.log[0], level="INFO")  # pyright: ignore
logger.add(sys.stderr, level="WARNING")  # Also log warnings to stderr

with logger.catch():
    logger.info("Starting MAGeCK design matrix generation")

    # Get parameters from snakemake
    metadata_file = snakemake.input.metadata  # pyright: ignore
    output_file = snakemake.output.design_matrix  # pyright: ignore
    design_column = snakemake.params.get("design_column", "deseq2")  # pyright: ignore

    logger.info(f"Reading metadata from: {metadata_file}")

    # Read metadata
    df = pd.read_csv(metadata_file, sep="\t")

    # Check if the design column exists
    if design_column not in df.columns:
        logger.error(f"Design column '{design_column}' not found in metadata")
        available_cols = ", ".join(df.columns)
        logger.error(f"Available columns: {available_cols}")
        raise ValueError(
            f"Design column '{design_column}' not found in metadata. "
            f"Available columns: {available_cols}"
        )

    # Check for missing values
    if df[design_column].isna().any():
        n_missing = df[design_column].isna().sum()
        logger.warning(
            f"Found {n_missing} samples with missing values in '{design_column}' column"
        )
        logger.warning("Samples with missing values will use 'unknown' as their group")
        df[design_column] = df[design_column].fillna("unknown")

    # Log the experimental groups found
    groups = df[design_column].unique()
    logger.info(f"Found {len(groups)} experimental groups: {', '.join(map(str, groups))}")

    # Count samples per group
    group_counts = df[design_column].value_counts()
    logger.info("Sample distribution across groups:")
    for group, count in group_counts.items():
        logger.info(f"  {group}: {count} samples")

    # Create design matrix using one-hot encoding
    logger.info("Creating design matrix with one-hot encoding")
    design_df = pd.get_dummies(df[design_column], prefix="", prefix_sep="")

    # Add baseline (intercept) column as the first column
    design_df.insert(0, "baseline", 1)

    # Use sample names as index (first column of metadata, typically 'samplename' or 'sample')
    sample_col = df.columns[0]
    design_df.index = df[sample_col]

    logger.info(f"Design matrix shape: {design_df.shape} (samples x conditions)")
    logger.info(f"Design matrix columns: {', '.join(design_df.columns)}")

    # Save design matrix
    logger.info(f"Saving design matrix to: {output_file}")
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    design_df.to_csv(output_file, sep="\t", index=True)

    # Log a preview of the design matrix
    logger.info("Design matrix preview (first 5 rows):")
    logger.info(f"\n{design_df.head().to_string()}")

    logger.info("MAGeCK design matrix generation completed successfully")
