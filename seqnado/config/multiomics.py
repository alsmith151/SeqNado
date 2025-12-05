from pathlib import Path


def find_assay_configs(directory: Path) -> tuple[dict, dict]:
    """Validate the existence of the config and metadata files for each assay."""

    config_files = {}
    metadata_files = {}

    for config_file in Path(directory).glob("config_*.yaml"):
        assay_name = config_file.stem.replace("config_", "")
        metadata_file = Path(directory) / f"metadata_{assay_name}.csv"

        if not metadata_file.exists():
            raise FileNotFoundError(
                f"Missing metadata file: {metadata_file}\n"
                f"Please create it using 'seqnado design {assay_name}'"
            )

        config_files[assay_name] = {"path": str(config_file)}
        metadata_files[assay_name] = str(metadata_file)

    return config_files, metadata_files
