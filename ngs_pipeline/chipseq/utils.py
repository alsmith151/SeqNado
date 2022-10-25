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


def get_sample_dataframe(files: List) -> pd.DataFrame:

    df = pd.DataFrame(files, columns=["fn"])

    df[["sample", "read"]] = (
        df["fn"].apply(str).str.extract("(?!.*/)?(.*)_.*_R?([12]).fastq.gz")
    )

    df["antibody"] = df["fn"].astype(str).str.split("_").str[-2]

    df = (
        df.pivot(columns="read", index=["sample", "antibody"])
        .droplevel(level=0, axis=1)
        .reset_index()
    )

    df_input = df.loc[df["antibody"].str.lower().str.contains("input")]
    df_input = df_input.assign(control=df_input["sample"] + "_" + df_input["antibody"])

    df_ip = df.loc[~df["antibody"].str.lower().str.contains("input")]

    df = df_ip.merge(df_input[["sample", "control"]], on="sample")

    return df


def get_pipeline_tools(config: Dict) -> Dict:

    return dict(
        homer="homer" in config["pileup_method"]
        or "homer" in config["peak_calling_method"],
        deeptools="deeptools" in config["pileup_method"],
        macs="macs" in config["peak_calling_method"],
        lanceotron="lanceotron" in config["peak_calling_method"],
    )
