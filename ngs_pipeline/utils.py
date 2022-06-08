
from typing import Dict, List
import os
import pathlib
import pandas as pd
import numpy as  np


def is_on(param: str) -> bool:
    """
    Returns True if parameter in "on" values
    On values:
        - true
        - t
        - on
        - yes
        - y
        - 1
    """
    values = ["true", "t", "on", "yes", "y", "1"]
    if str(param).lower() in values:
        return True
    else:
        return False



def is_off(param: str):
    """Returns True if parameter in "off" values"""
    values = ["", "None", "none", "F", "f"]
    if str(param).lower() in values:
        return True
    else:
        return False


def is_none(param: str) -> bool:
    """Returns True if parameter is none"""
    values = ["", "none"]
    if str(param).lower() in values:
        return True
    else:
        return False


def convert_empty_yaml_entry_to_string(param: str) -> str:
    """
    Converts empty yaml entries to string
    """
    if is_none(param):
        return ""
    else:
        return param

def format_config_dict(config: Dict) -> Dict:
    """
    Formats the config dictionary to ensure that all entries are strings.

    """
    for key, value in config.items():
        config[key] = convert_empty_yaml_entry_to_string(value)

    return config


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


def get_fastq_files(path: str, recursive=False) -> pd.DataFrame:

    data = pathlib.Path(path).glob("**/*.fastq.gz") if recursive else pathlib.Path(path).glob("*.fastq.gz")
    df = pd.DataFrame(data=data, columns=["fn"])
    df = df.assign(basename=lambda df: df["fn"].apply(lambda p: p.name))

    df["paired_or_single"] = df["basename"].str.match(r"(.*)_R?[12].fastq(.gz)?")
    df["paired_or_single"] = np.where(
        df["paired_or_single"] == True, "paired", "single"
    )
    return df