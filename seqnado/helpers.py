from typing import Dict
import pathlib
import numpy as np


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
