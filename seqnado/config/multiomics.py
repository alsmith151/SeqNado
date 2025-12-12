from pathlib import Path
from typing import Annotated

from pydantic import BaseModel, BeforeValidator, field_validator


def none_str_to_none(v):
    """Convert string 'none' to None."""
    if isinstance(v, str) and v.strip().lower() == "none":
        return None
    return v


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
    output_dir: Annotated[
        str | None,
        BeforeValidator(none_str_to_none)
    ] = "seqnado_output/"

    @property
    def summary_report(self) -> str:
        """Path to the multi-assay summary report."""
        return str(Path(self.output_dir) / "multi_assay_summary.txt")

    @property
    def heatmap(self) -> str:
        """Path to the multi-assay heatmap PDF."""
        return str(Path(self.output_dir) / "multiomics" / "heatmap" / "heatmap.pdf")

    @property
    def metaplot(self) -> str:
        """Path to the multi-assay metaplot PDF."""
        return str(Path(self.output_dir) / "multiomics" / "heatmap" / "metaplot.pdf")

    @property
    def all_outputs(self) -> list[str]:
        """Get all multiomics output files."""
        return [
            self.summary_report,
            self.heatmap,
            self.metaplot,
        ]
