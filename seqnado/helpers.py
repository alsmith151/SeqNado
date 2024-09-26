from typing import Dict, Union, Optional, List, Tuple, Any
import pathlib
import numpy as np
import shlex

from loguru import logger


from seqnado.design import Design, DesignIP


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
        options = [o for i, o in enumerate(options) if i not in [cores_flag, cores_flag + 1]]
    except ValueError:
        try:
            cores_flag = options.index("--cores")
            cores = int(options[cores_flag + 1])
            options = [o for i, o in enumerate(options) if i not in [cores_flag, cores_flag + 1]]
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
    from loguru import logger

    try:
        apptainer_flag = options.index("--apptainer-args")
        apptainer_args = options[apptainer_flag + 1]
        options = [o for i, o in enumerate(options) if i not in [apptainer_flag, apptainer_flag + 1]]
    except ValueError:
        apptainer_args = ""

    return options, apptainer_args


def define_memory_requested(attempts: int = 1, initial_value: int  = 1, scale: float = 1) -> str:
    """
    Define the memory requested for the job.
    """
    memory = int(initial_value) * 2 ** (int(attempts) - 1)
    memory = memory * float(scale)
    return f"{memory}G"

def define_time_requested(attempts: int = 1, initial_value: int = 1, scale: float = 1) -> str:
    """
    Define the time requested for the job.

    Base time is 1 hour.
    """
    time = int(initial_value) * 2 ** (int(attempts) - 1)
    time = time * float(scale)
    return f"{time}h"


def symlink_file(
    output_dir: pathlib.Path, source_path: pathlib.Path, new_file_name: str
):
    """
    Create a symlink in the output directory with the new file name.
    """

    new_path = output_dir / new_file_name
    if not new_path.exists() and source_path.is_file():
        # logger.debug(f"Symlinking {source_path} to {output_dir / new_file_name}")
        if str(source_path) in [".", "..", "", None, "None"]:
            logger.warning(f"Source path is empty for {new_file_name}. Will not symlink.")

        else:
            new_path.symlink_to(source_path.resolve())
            # logger.debug(f"Symlinked {source_path} to {output_dir / new_file_name} successfully.")


def symlink_fastq_files(
    design: Union[Design, DesignIP], output_dir: str = "seqnado_output/fastqs/"
):
    """
    Symlink the fastq files to the output directory.
    """
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(design, Design):
        for fastq_set in design.fastq_sets:
            if fastq_set.is_paired:
                symlink_file(output_dir, fastq_set.r1.path, f"{fastq_set.name}_1.fastq.gz")
                symlink_file(output_dir, fastq_set.r2.path, f"{fastq_set.name}_2.fastq.gz")
            else:
                symlink_file(output_dir, fastq_set.r1.path, f"{fastq_set.name}.fastq.gz")

    elif isinstance(design, DesignIP):
        for experiment in design.experiments:

            if experiment.fastqs_are_paired:
                symlink_file(
                    output_dir, experiment.ip.r1.path, f"{experiment.ip.name}_{experiment.ip_performed}_1.fastq.gz"
                )
                symlink_file(
                    output_dir, experiment.ip.r2.path, f"{experiment.ip.name}_{experiment.ip_performed}_2.fastq.gz"
                )
            else:
                symlink_file(output_dir, experiment.ip.r1.path, f"{experiment.ip.name}_{experiment.ip_performed}.fastq.gz")

            # Control files
            if experiment.has_control:
                if experiment.control.is_paired:
                    symlink_file(
                        output_dir,
                        experiment.control.r1.path,
                        f"{experiment.control.r1.sample_base_without_ip}_{experiment.control_performed}_1.fastq.gz",
                    )
                    symlink_file(
                        output_dir,
                        experiment.control.r2.path,
                        f"{experiment.control.r1.sample_base_without_ip}_{experiment.control_performed}_2.fastq.gz",
                    )
                else:
                    symlink_file(
                        output_dir,
                        experiment.control.r1.path,
                        f"{experiment.control.r1.sample_base_without_ip}_{experiment.control_performed}.fastq.gz",
                    )


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
    values = ["", "none", "f", "n", "no", "false", "0"]
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
        if isinstance(value, dict):
            config[key] = format_config_dict(value)
        else:
            entry = convert_empty_yaml_entry_to_string(value)

            if is_on(entry):
                config[key] = True
            elif is_off(entry):
                config[key] = False
            elif is_none(entry):
                config[key] = False
            else:
                config[key] = entry

    return config


def has_bowtie2_index(prefix: str) -> bool:
    """
    Checks if bowtie2 index is present.
    """

    path_prefix = pathlib.Path(prefix).resolve()
    path_dir = path_prefix.parent
    path_prefix_stem = path_prefix.stem

    bowtie2_indices = list(path_dir.glob(f"{path_prefix_stem}*.bt2"))

    if len(bowtie2_indices) > 0:
        return True


def check_options(value: object):
    if value in [None, np.nan, ""]:
        return ""
    elif is_off(value):
        return ""
    else:
        return value

def pepe_silvia():
    print("PEPE SILVIA")
    _pepe_silvia = (
        "https://relix.com/wp-content/uploads/2017/03/tumblr_o16n2kBlpX1ta3qyvo1_1280.jpg"
    )
    return _pepe_silvia


def get_group_for_sample(wildcards, design: Union[Design, DesignIP], strip: str = ""):
    from seqnado.design import NormGroups

    norm_groups = NormGroups.from_design(design, include_controls=True)

    try:
        group = norm_groups.get_sample_group(wildcards.sample.strip(strip))
        return group
    except KeyError:
        # logger.error(f"Sample {wildcards.sample} not found in normalisation groups.")
        raise KeyError(f"Sample {wildcards.sample} not found in normalisation groups.")

def get_scale_method(config: Dict) -> Optional[str]:
    """
    Returns the scale method based on the config.
    """

    if config["spikein"]:
       method = "spikein"
    elif config["scale"]:
        method = "csaw"
    else:
        method = None
    
    return method


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


def get_genomes(): 

    import json
    import sys

    seqnado_config_dir = pathlib.Path('~/.config/seqnado').expanduser()
    genome_config = seqnado_config_dir / 'genomes.json'

    if genome_config.exists():
        logger.debug(f"Genome config found at {genome_config}")
        with open(genome_config, 'r') as f:
            genome_values = json.load(f)
        
        if genome_values.get('dm6'):
            if "PATH" in  genome_values['dm6']['bt2_indices']:
                logger.error(f'Template genome file found. Please update the genome file {genome_config} with the correct paths.')
                sys.exit(1)
    
    return genome_values



        






    # with open(os.path.join(template_dir, "preset_genomes.json"), "r") as f:
    #     genome_values = json.load(f)
