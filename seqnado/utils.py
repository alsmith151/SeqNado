from typing import Dict, List, Union
import os
import pathlib
import pandas as pd
import numpy as np
import snakemake
import hashlib
from collections import defaultdict


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


def get_fastq_files(path: str, recursive=False) -> pd.DataFrame:

    files = (
        pathlib.Path(path).glob("**/*.fastq.gz")
        if recursive
        else pathlib.Path(path).glob("*.fastq.gz")
    )
    return files


def get_singularity_command(workflow: snakemake.Workflow, command: str):
    """
    Runs a command in a singularity container.
    """

    container_url = workflow.global_container_img
    container_dir = workflow.persistence.container_img_path

    md5 = hashlib.md5()
    md5.update(container_url.encode())
    container_hash = md5.hexdigest()

    container = os.path.join(container_dir, container_hash + ".simg")

    command = command.replace(r"\n", r" \\")
    command = f"bash -c '{command}'"

    return f"singularity exec -H $PWD {workflow.singularity_args} {container} {command}"


class GenericFastqSamples:
    def __init__(self, design):

        # Expected columns: sample, fq1, fq2
        self.design = design
        self.design = self.design.assign(
            paired=(~self.design[["fq1", "fq2"]].isna().any(axis=1))
        )

    @classmethod
    def from_files(cls, files: List[Union[pathlib.Path, str]]) -> "GenericFastqSamples":

        df = pd.DataFrame(files, columns=["fn"])

        df[["sample", "read"]] = (
            df["fn"].apply(str).str.extract("(?!.*/)?(.*).*_R?([12]).fastq.gz")
        )

        df["sample"] = df["sample"].apply(
            lambda p: pathlib.Path(p).name
            if isinstance(p, pathlib.Path)
            else os.path.basename(p)
        )
        df["read"] = "fq" + df["read"]

        df = (
            df.pivot(columns="read", index=["sample"])
            .droplevel(level=0, axis=1)
            .reset_index()
        )

        return cls(design=df)

    @property
    def fastq_files(self):
        return sorted([*self.design["fq1"], *self.design["fq2"]])

    @property
    def sample_names_all(self):
        return self.design["sample"].to_list()

    @property
    def translation(self):
        fq_translation = {}
        for sample in self.design.itertuples():
            for read in [1, 2]:
                fq_translation[f"{sample.sample}_{read}.fastq.gz"] = os.path.realpath(
                    str(getattr(sample, f"fq{read}"))
                )

        return fq_translation


def check_options(value: object):
    if value in [None, np.nan, ""]:
        return ""
    else:
        return value


def translate_fq_files(wc, samples: GenericFastqSamples, paired: bool=False):

    if paired:
        return {"fq1": samples.translation[f"{wc.sample}_1.fastq.gz"],
                "fq2": samples.translation[f"{wc.sample}_2.fastq.gz"]}
    else:
        return {"fq": samples.translation[f"{wc.sample}_{wc.read}.fastq.gz"]}

def get_fq_filestem(wc, samples: GenericFastqSamples):
    fn = samples.translation[f"{wc.sample}_{wc.read}.fastq.gz"]
    basename = os.path.basename(fn)
    return os.path.splitext(basename.replace(".gz", ""))[0]


def pair_treatment_and_control_for_peak_calling(wc, samples, assay, filetype):

    if assay == "ChIP":

        df_design_sample = samples.loc[(samples["sample"] == wc.sample) & (samples["antibody"] == wc.antibody)]
        if df_design_sample.empty:
            raise Exception(f"Could not find sample {wc.sample} with antibody {wc.antibody} in design file")

        filetype_to_dir_mapping = {"tag": "tag_dirs", "bigwig": "bigwigs/deeptools", "bam": "aligned"}
        filetype_to_extension_mapping = {"tag": "/", "bigwig": ".bigWig", "bam": ".bam"}
        
        extension_for_filetype = filetype_to_extension_mapping[filetype]
        directory_for_filetype = filetype_to_dir_mapping[filetype]
        
        treatment = f"seqnado_output/{directory_for_filetype}/{wc.sample}_{wc.antibody}{extension_for_filetype}"
        control = f"seqnado_output/{directory_for_filetype}/{df_design_sample.iloc[0]['control']}{extension_for_filetype}"


        files =  {"treatment": treatment, "control": control}

    else:
        files = {"treatment": treatment, "control": ""}
    
    return files

