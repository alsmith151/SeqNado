import json
import logging
import pathlib
import sys
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from loguru import logger

from seqnado.inputs import SampleCollection, SampleCollectionForIP, ScaleMethod

FILETYPE_TO_DIR_MAPPING = {
    "tag": "tag_dirs",
    "bigwig": "bigwigs/deeptools",
    "bam": "aligned",
}

FILETYPE_TO_EXTENSION_MAPPING = {"tag": "/", "bigwig": ".bigWig", "bam": ".bam"}


def extract_cores_from_options(options: List[str]) -> Tuple[List[str], int]:
    """
    Extract the number of cores from the snakemake options.
    """

    try:
        cores_flag = options.index("-c")
        cores = int(options[cores_flag + 1])
        options = [
            o for i, o in enumerate(options) if i not in [cores_flag, cores_flag + 1]
        ]
    except ValueError:
        try:
            cores_flag = options.index("--cores")
            cores = int(options[cores_flag + 1])
            options = [
                o
                for i, o in enumerate(options)
                if i not in [cores_flag, cores_flag + 1]
            ]
        except ValueError:
            cores = 1
            logger.warning("No core flag provided. Defaulting to 1 core.")
    except IndexError:
        cores = 1
        options = [o for i, o in enumerate(options) if i not in [cores_flag]]
        logger.warning("Core flag provided but no value given. Defaulting to 1 core.")

    return options, cores


def extract_apptainer_args(options: List[str]) -> Tuple[List[str], str]:
    """
    Extract the apptainer arguments from the snakemake options.
    """
    try:
        apptainer_flag = options.index("--apptainer-args")
        apptainer_args = options[apptainer_flag + 1]
        options = [
            o
            for i, o in enumerate(options)
            if i not in [apptainer_flag, apptainer_flag + 1]
        ]
    except ValueError:
        apptainer_args = ""

    return options, apptainer_args

def define_memory_requested(
    attempts: int = 1,
    initial_value: int = 1,
    scale: float = 1
) -> str:
    """
    Define the memory requested for the job, 
    returns a string like "4G" avoiding decimals for qualimap.
    """
    mem_value = int(initial_value) * 2 ** (int(attempts) - 1)
    mem_value = int(mem_value * float(scale))
    return f"{mem_value}G"



def define_time_requested(
    attempts: int = 1, initial_value: int = 1, scale: float = 1
) -> str:
    """
    Define the time requested for the job.

    Base time is 1 hour.
    """
    time = int(initial_value) * 2 ** (int(attempts) - 1)
    time = time * float(scale)
    return f"{time}h"


def pepe_silvia():
    print("PEPE SILVIA")
    _pepe_silvia = "https://relix.com/wp-content/uploads/2017/03/tumblr_o16n2kBlpX1ta3qyvo1_1280.jpg"
    return _pepe_silvia


def get_group_for_sample(wildcards, design: Union[SampleCollection, SampleCollectionForIP], strip: str = ""):
    from seqnado.inputs import SampleGroups

    norm_groups = SampleGroups.from_sample_collection(design, include_controls=True)

    try:
        group = norm_groups.get_sample_group(wildcards.sample.strip(strip))
        return group
    except KeyError:
        # logger.error(f"Sample {wildcards.sample} not found in normalisation groups.")
        raise KeyError(f"Sample {wildcards.sample} not found in normalisation groups.")


def get_scale_method(config: Dict) -> List[str]:
    """
    Returns the scale method based on the config.
    """

    method = [ScaleMethod.unscaled]

    if config.get("spikein"):
        method.append(ScaleMethod.spikein)
    elif config.get("scale"):
        method.append(ScaleMethod.csaw)
    return [m.value for m in method]

def remove_unwanted_run_files():
    import glob
    import os
    import shutil

    slurm_files = glob.glob("slurm-*.out")
    sps_files = glob.glob("sps-*")
    simg_files = glob.glob("*.simg")

    for fn in [*slurm_files, *sps_files, *simg_files]:
        try:
            if not os.path.isdir(fn):
                os.remove(fn)
            else:
                shutil.rmtree(fn)

        except Exception as e:
            print(e)


def run_batch_job_on_error(email: str):
    """
    Run a batch job with slurm to send an email on error.
    """

    import subprocess

    slurm_script = f"""#!/bin/bash
    #SBATCH --job-name=seqnado_error_notification
    #SBATCH --mail-type=END
    #SBATCH --mail-user={email}

    echo "An error occurred in the job. Please check the logs for more details."
    """

    with open("error_notification.sh", "w") as f:
        f.write(slurm_script)

    try:
        subprocess.run(["sbatch", "error_notification.sh"], check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to submit slurm job: {e}")


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
        df.columns = bed_columns[:len(df.columns)]
        # Ensure we have at least the minimum required columns
        if "Name" not in df.columns:
            df["Name"] = df["Chromosome"] + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
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
