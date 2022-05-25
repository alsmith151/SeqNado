import os
from typing import Dict, List
import pathlib
import pybedtools
import pandas as pd
import re
import glob
import numpy as np


def set_up_chromsizes(config: Dict):
    """
    Ensures that genome chromsizes are present.

    If chromsizes are not provided this function attempts to download them from UCSC.
    The P.PARAMS dictionary is updated with the location of the chromsizes.

    """

    try:
        config["genome"]["name"]
    except KeyError:
        raise "Genome name has not been provided."

    if config["genome"].get("chrom_sizes"):
        pass

    elif os.path.exists("chrom_sizes.txt.tmp"):
        config["genome"]["chrom_sizes"] = "chrom_sizes.txt.tmp"

    else:
        from pybedtools.helpers import get_chromsizes_from_ucsc

        get_chromsizes_from_ucsc(config["genome"]["name"], "chrom_sizes.txt.tmp")
        config["genome"]["chrom_sizes"] = "chrom_sizes.txt.tmp"

    print(f"Obtained chromosome sizes for genome {config['genome']['name']}")


def symlink_fastq_files(sample_info: pd.DataFrame):

    try:
        os.mkdir("fastq")
    except FileExistsError:
        pass

    for fq in sample_info.itertuples():
        full_path = pathlib.Path(fq.fn).absolute()
        symlink_path = f"fastq/{fq.basename}"

        if not os.path.exists(symlink_path):
            os.symlink(full_path, symlink_path)


def get_fastq_files(path: str) -> pd.DataFrame:
    df = pd.DataFrame(data=pathlib.Path(".").glob("*.fastq.gz"), columns=["fn"])
    df = df.assign(basename=lambda df: df["fn"].apply(lambda p: p.name))

    df["paired_or_single"] = df["basename"].str.match(r"(.*)_R?[12].fastq(.gz)?")
    df["paired_or_single"] = np.where(
        df["paired_or_single"] == True, "paired", "single"
    )
    return df


def sample_names_follow_convention(
    df: pd.DataFrame, name_column: str = "basename"
) -> bool:
    naming_pattern_paired = r"(.*)_(.*)_R?[12].fastq(.gz)?"
    naming_pattern_single = r"(.*)_(.*).fastq(.gz)?"

    return (
        df[name_column].str.match(naming_pattern_paired)
        | df[name_column].str.match(naming_pattern_single)
    ).all()


def get_sample_attributes(df: pd.DataFrame):

    if not df.empty:
        if not sample_names_follow_convention(df):
            raise ValueError(
                "Sample names do not follow the naming convention. Use a design matrix to specify inputs"
            )

    df_paired = df.query("paired_or_single == 'paired'")
    df_single = df.query("paired_or_single == 'single'")

    df_paired = df_paired.join(
        df_paired["basename"].str.extract(
            r"(?P<sample_name>.*?)_(?P<antibody_or_input>.*?)_(?P<read_number>[12]).fastq.gz"
        )
    )

    df_single = df_single.join(
        df_single["basename"]
        .str.extract(r"(?P<sample_name>.*?)_(?P<antibody_or_input>.*?).fastq.gz")
        .assign(read_number=0)
    )

    return pd.concat([df_paired, df_single]).assign(
        fn=df["fn"].astype(str),
        basename=df["basename"],
        paired=df["paired_or_single"].str.contains("paired"),
    )


def has_ambigous_inputs(df: pd.DataFrame) -> bool:
    return (
        df.drop_duplicates(subset=["sample_name", "antibody", "input"])
        .duplicated(subset=["sample_name", "input"])
        .any()
    )


def pair_inputs_with_samples(df: pd.DataFrame):

    df_input = df.query("antibody_or_input.str.contains('input', case=False)")
    df_ip = df.query("~antibody_or_input.str.contains('input', case=False)")

    df_ip_with_input = df_ip.merge(
        df_input[["sample_name", "antibody_or_input"]], on="sample_name", how="left"
    ).rename(
        columns={"antibody_or_input_x": "antibody", "antibody_or_input_y": "input"}
    )

    if not has_ambigous_inputs(df_ip_with_input):
        return df_ip_with_input
    else:
        raise ValueError(
            "Samples are paired with more than one input. Use a design matrix to specify inputs"
        )
