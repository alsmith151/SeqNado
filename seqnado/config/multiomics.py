from pathlib import Path
from typing import Annotated

from pydantic import BaseModel, BeforeValidator


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


class MultiomicsOutput(BaseModel):
    output_dir: Annotated[str | None, BeforeValidator(none_str_to_none)] = (
        "seqnado_output/"
    )

    @property
    def summary_report(self) -> str:
        """Path to the multiomics summary report."""
        return str(Path(self.output_dir) / "multiomics_summary.txt")

    @property
    def heatmap(self) -> str:
        """Path to the multiomics heatmap PDF."""
        return str(Path(self.output_dir) / "multiomics" / "heatmap" / "heatmap.pdf")

    @property
    def metaplot(self) -> str:
        """Path to the multiomics metaplot PDF."""
        return str(Path(self.output_dir) / "multiomics" / "heatmap" / "metaplot.pdf")
    
    @property
    def dataset(self) -> str:
        """Get the output directory."""
        return str(Path(self.output_dir) / "multiomics" / "dataset" / "dataset_regions.h5ad")

    @property
    def all_outputs(self) -> list[str]:
        """Get all multiomics output files."""
        return [
            self.summary_report,
            self.heatmap,
            self.metaplot,
            self.dataset,
        ]
