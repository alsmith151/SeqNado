from pathlib import Path
from typing import Annotated

from pydantic import BaseModel, BeforeValidator, Field
from seqnado import Assay


def none_str_to_none(v):
    """Convert string 'none' to None."""
    if isinstance(v, str) and v.strip().lower() == "none":
        return None
    return v


def get_assay_bigwigs(wildcards, rules, ASSAYS) -> list[str]:
    """Get all bigwigs from assay-specific 'all' rules."""
    bigwigs = []
    for assay in ASSAYS:
        rule_name = f"{assay}_all"
        inputs = getattr(rules, rule_name).input

        # Handle both single files and collections of files
        if isinstance(inputs, str):
            if inputs.endswith(".bigWig"):
                bigwigs.append(inputs)
        else:
            # InputFiles is iterable
            for file in inputs:
                if isinstance(file, str) and file.endswith(".bigWig"):
                    bigwigs.append(file)
    return bigwigs


def find_assay_configs(directory: Path) -> tuple[dict, dict]:
    """Validate the existence of the config and metadata files for each assay.

    Note: Skips config_multiomics.yaml as it's not an assay-specific config.
    """

    config_files = {}
    metadata_files = {}

    for config_file in Path(directory).glob("config_*.yaml"):
        assay_name = config_file.stem.replace("config_", "")

        # Skip multiomics config - it's not an assay config
        if assay_name == "multiomics":
            continue

        metadata_file = Path(directory) / f"metadata_{assay_name}.csv"

        if not metadata_file.exists():
            raise FileNotFoundError(
                f"Missing metadata file: {metadata_file}\n"
                f"Please create it using 'seqnado design {assay_name}'"
            )

        config_files[assay_name] = {"path": str(config_file)}
        metadata_files[assay_name] = str(metadata_file)

    return config_files, metadata_files