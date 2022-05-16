import os
from typing import Dict, List
import pathlib
import pybedtools
import pandas as pd
import re
import glob


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


def get_fastq_files(path: os.PathLike) -> pd.DataFrame:
    fq_list = glob.glob(path)
    sample_info = pd.DataFrame(fq_list, columns=["fq_files"])
    sample_info["fq"] = [os.path.basename(fn) for fn in sample_info["fq_files"]]

    sample_info["is_single_end"] = sample_info["fq"].str.match(r"(.*?).fastq(.gz)?")
    sample_info["is_paired_end"] = sample_info["fq"].str.match(
        r"(.*)_R?[12].fastq(.gz)?"
    )

    sample_info_single = sample_info.loc[sample_info["is_single_end"]]
    sample_info_paired = sample_info.loc[sample_info["is_paired_end"]]

    sample_info_single["sample_name"] = sample_info_single["fq"].str.extract(
        "(.*?).fastq(?!.gz)?"
    )
    sample_info_single["read"] = 0

    sample_info_paired["sample_name"] = sample_info_paired["fq"].str.extract(
        "(.*?)_R?\d.fastq(?!.gz)?"
    )
    sample_info_paired["read"] = (
        sample_info_paired["fq"].str.extract(".*?_R?(\d).fastq(?!.gz)?").astype(int)
    )

    return pd.concat([sample_info_single, sample_info_paired])


def extract_chip_sample_info(sample_info: pd.DataFrame):

    # Ensure that all samples match the naming convention
    if (
        not sample_info["fq"]
        .str.match("(?P<sample>.*)_(?P<antibody>.*?)(_R?[12])?.fastq(.gz)?")
        .all()
    ):
        raise ValueError("Sample names do not match the naming convention.")

    # Extract sample names and antibody names
    sample_info[["sample_name", "antibody"]] = sample_info["fq"].str.extract(
        "(?P<sample>.*?)_(?P<antibody>.*?)(?!_R?[12])?.fastq(?!.gz)?"
    )
    sample_info["antibody"] = sample_info["antibody"].str.split("_[12]").str[0]

    sample_info["is_input"] = sample_info["sample_name"].str.contains(
        "input", case=False
    )

    return {
        "info": sample_info,
        "ip": sample_info.query("is_input == False")["fq_files"].to_list(),
        "input": sample_info.query("is_input == True")["fq_files"].to_list(),
    }


def symlink_fastq_files(sample_info: pd.DataFrame):

    try:
        os.mkdir("fastq")
    except FileExistsError:
        pass

    for fq in sample_info.itertuples():
        full_path = pathlib.Path(fq.fq_files).absolute()
        symlink_path = f"fastq/{fq.sample_name}_{fq.antibody}_R{fq.read}.fastq.gz"

        if not os.path.exists(symlink_path):
            os.symlink(full_path, symlink_path)
