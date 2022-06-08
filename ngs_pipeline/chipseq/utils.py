import os
from typing import Dict, List
import pathlib
import pybedtools
import pandas as pd
import re
import glob
import numpy as np
from typing import List, Literal

pd.set_option("mode.chained_assignment", None)


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
        .dropna(subset=["input"])
        .duplicated(subset=["ip", "input"])
        .any()
    )


def pair_inputs_with_samples(df: pd.DataFrame):

    df_input = df.query("antibody_or_input.str.contains('input', case=False)")
    df_ip = df.query("~antibody_or_input.str.contains('input', case=False)")

    df_ip["ip"] = df_ip["sample_name"] + "_" + df_ip["antibody_or_input"]
    df_input["input"] = df_input["sample_name"] + "_" + df_input["antibody_or_input"]

    df_ip_with_input = (
        df_ip.merge(df_input[["sample_name", "input"]], on="sample_name", how="left")
        .rename(columns={"antibody_or_input": "antibody"})
        .drop_duplicates(["basename"])
        .sort_values("basename")
        .reset_index(drop=True)
    )

    if not has_ambigous_inputs(df_ip_with_input):
        return df_ip_with_input
    else:
        raise ValueError(
            "Samples are paired with more than one input. Use a design matrix to specify inputs"
        )

    return df_ip_with_input


def get_pipeline_tools(config: Dict) -> Dict:

    return dict(
        homer="homer" in config["pileup_method"]
        or "homer" in config["peak_calling_method"],
        deeptools="deeptools" in config["pileup_method"],
        macs="macs" in config["peak_calling_method"],
        lanceotron="lanceotron" in config["peak_calling_method"],
    )
