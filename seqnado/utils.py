from typing import Dict, List, Union, Tuple, Literal
import os
import pathlib
import pandas as pd
import numpy as np
import snakemake
import hashlib
from collections import defaultdict

from snakemake.io import expand


FILETYPE_TO_DIR_MAPPING = {
            "tag": "tag_dirs",
            "bigwig": "bigwigs/deeptools",
            "bam": "aligned",
        }

FILETYPE_TO_EXTENSION_MAPPING = {"tag": "/", "bigwig": ".bigWig", "bam": ".bam"}




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
    values = ["", "None", "none", "F", "f", "n", "no"]
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

        if isinstance(config[key], str):
            if is_on(config[key]):
                config[key] = True
            elif is_off(config[key]):
                config[key] = False
            elif is_none(config[key]):
                config[key] = None

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


def translate_fq_files(wc, samples: GenericFastqSamples, paired: bool = False):

    if paired:
        return {
            "fq1": samples.translation[f"{wc.sample}_1.fastq.gz"],
            "fq2": samples.translation[f"{wc.sample}_2.fastq.gz"],
        }
    else:
        return {"fq": samples.translation[f"{wc.sample}_{wc.read}.fastq.gz"]}


def get_fq_filestem(wc, samples: GenericFastqSamples):
    fn = samples.translation[f"{wc.sample}_{wc.read}.fastq.gz"]
    basename = os.path.basename(fn)
    return os.path.splitext(basename.replace(".gz", ""))[0]



def get_treatment_file(wc, assay, filetype):

    extension_for_filetype = FILETYPE_TO_EXTENSION_MAPPING[filetype]
    directory_for_filetype = FILETYPE_TO_DIR_MAPPING[filetype]

    treatment = f"seqnado_output/{directory_for_filetype}/{wc.treatment}{extension_for_filetype}"
    
    return treatment

def get_control_file(wc, design,  assay, filetype):

    extension_for_filetype = FILETYPE_TO_EXTENSION_MAPPING[filetype]
    directory_for_filetype = FILETYPE_TO_DIR_MAPPING[filetype]

    if assay == "ChIP":

        df = design.assign(treatment = lambda df: df[["sample", "antibody"]]["sample"].str.cat(df["antibody"], sep="_"))
        sample_row = df.query("treatment == @wc.treatment").iloc[0]
        has_control = not pd.isna(sample_row["control"])

        if has_control:
            control = f"seqnado_output/{directory_for_filetype}/{sample_row['control']}{extension_for_filetype}"
        else:
            control = "NA"
    
    return control


def define_output_files(
    assay: Literal["ChIP", "ATAC", "RNA", "SNP"],
    sample_names: list = None,
    pileup_method: list = None,
    peak_calling_method: list = None,
    make_bigwigs: bool = False,
    call_peaks: bool = False,
    make_heatmaps: bool = False,
    make_ucsc_hub: bool = False,
    call_snps: bool = False,
    annotate_snps: bool = False,
    **kwargs,
) -> list:
    """Define output files for the pipeline"""

    analysis_output = ["seqnado_output/qc/full_qc_report.html"]
    assay_output = []

    if make_ucsc_hub:
        hub_dir = kwargs["ucsc_hub_details"].get("directory")
        hub_name = kwargs["ucsc_hub_details"].get("name")
        hub_file = os.path.join(hub_dir,  f"{hub_name}.hub.txt")
        analysis_output.append(hub_file)


    if assay in ["ChIP", "ATAC"]:

        if make_bigwigs and pileup_method:
            assay_output.extend(
                expand(
                    "seqnado_output/bigwigs/{method}/{sample}.bigWig",
                    sample=sample_names,
                    method=pileup_method,
                )
            )

        if call_peaks and peak_calling_method:
            if assay == "ChIP":
                assay_output.extend(
                    expand(
                        "seqnado_output/peaks/{method}/{ip}.bed",
                        ip=kwargs["sample_names_ip"],
                        method=peak_calling_method,
                    )
                )
            else:
                assay_output.extend(
                    expand(
                        "seqnado_output/peaks/{method}/{sample}.bed",
                        sample=sample_names,
                        method=peak_calling_method,
                    )
                )

        if make_heatmaps:
            assay_output.extend(
                expand(
                    "seqnado_output/heatmap/{method}/{sample}.png",
                    sample=sample_names,
                    method=pileup_method,
                )
            )

    elif assay == "RNA":

        if make_bigwigs and pileup_method:
            assay_output.extend(
                expand(
                    "seqnado_output/bigwigs/{method}/{sample}_{strand}.bigWig",
                    sample=sample_names,
                    method=pileup_method,
                    strand=["plus", "minus"],
                )
            )
        
        if kwargs["run_deseq2"]:
            project_id = kwargs["DESeq2"].get("project_id")
            assay_output.append(f"DESeq2_{project_id}.html") 


        assay_output.extend(
            [
                "seqnado_output/feature_counts/read_counts.tsv",
                *expand(
                    "seqnado_output/aligned/{sample}.bam",
                    sample=sample_names,
                ),
            ]
        )

    elif assay == "SNP":

        if call_snps:
            assay_output.expand(
                "seqnado_output/variant/{sample}_filtered.anno.vcf.gz",
                sample=sample_names,
            )

        if annotate_snps:
            assay_output.append(
                expand(
                    "seqnado_output/variant/{sample}_filtered.stats.txt",
                    sample=sample_names,
                ),
            )

    analysis_output.extend(assay_output)


    return analysis_output
