import json
import re
import tarfile
import time
from dataclasses import dataclass
from pathlib import Path

import pytest
import requests
import yaml


def load_yaml(path: Path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def load_json(path: Path):
    with open(path, "r") as f:
        return json.load(f)


def download_with_retry(url: str, dest: Path, max_retries: int = 3, timeout: int = 30):
    for attempt in range(max_retries):
        try:
            dest.parent.mkdir(
                parents=True, exist_ok=True
            )  # Ensure parent directory exists
            r = requests.get(url, stream=True, timeout=timeout)
            r.raise_for_status()
            with open(dest, "wb") as f:
                f.write(r.content)
            return
        except (requests.RequestException, OSError):
            if attempt < max_retries - 1:
                time.sleep(2**attempt)
                continue
            raise


def extract_tar(tar_path: Path, dest: Path, flatten: bool = True):
    with tarfile.open(tar_path) as tar:
        tar.extractall(path=dest)
    if flatten:
        nested = dest / dest.name
        if nested.exists() and nested.is_dir():
            for f in nested.iterdir():
                f.rename(dest / f.name)
            try:
                nested.rmdir()
            except Exception as e:
                print(f"[WARNING] Could not remove nested directory {nested}: {e}")


def setup_genome_config(
    genome_config_file: Path,
    star_index: Path,
    bt2_index: Path,
    chromsizes: Path,
    gtf: Path,
    blacklist: Path,
    genes_bed: Path,
    fasta: Path,
):
    """
    Write a nested genome_config.json for seqnado tests.
    """
    config_data = {
        "hg38": {
            "star_index": str(star_index) if star_index else None,
            "bt2_index": str(bt2_index) if bt2_index else None,
            "chromosome_sizes": str(chromsizes) if chromsizes else None,
            "gtf": str(gtf) if gtf else None,
            "blacklist": str(blacklist) if blacklist else None,
            "genes": str(genes_bed) if genes_bed else None,
            "fasta": str(fasta) if fasta else None,
        }
    }

    genome_config_file.parent.mkdir(parents=True, exist_ok=True)
    with open(genome_config_file, "w") as f:
        json.dump(config_data, f, indent=2)


# Fix path to point to tests/pipeline/config
def get_fastq_pattern(assay: str) -> str:
    """Get the FASTQ file pattern for a given assay."""
    patterns = {
        "atac": "atac_*.fastq.gz",
        "chip": "chip-rx_*.fastq.gz",
        "chip-rx": "chip-rx_*.fastq.gz",
        "rna": "rna_*.fastq.gz",
        "rna-rx": "rna-spikein-*.fastq.gz",
        "snp": "snp_*.fastq.gz",
        "cat": "chip-rx_*.fastq.gz",
        "meth": "meth-*.fastq.gz",
        "mcc": "mcc_*.fastq.gz",
    }
    pattern = patterns.get(assay)
    if not pattern:
        raise ValueError(f"Unsupported assay: {assay}")
    return pattern


# Utility function to parse assays from pytest config
def parse_assays(config: pytest.Config) -> list[str]:
    """Parse --assays option into list."""
    assays_str: str = config.getoption("--assays") or "chip"
    assays = [a.strip() for a in assays_str.split(",") if a.strip()]
    return assays or ["chip"]


# TestPaths dataclass for path management
@dataclass(frozen=True)
class TestPaths:
    """
    Container for all important test and project paths.
    All paths are created and ready for use.
    """

    repo: Path
    package: Path
    test_dir: Path
    test_data: Path
    workflow: Path
    test_profile: Path
    genome: Path
    fastq: Path


# Factory function for TestPaths
def make_test_paths(file: Path) -> TestPaths:
    # Search upwards for pyproject.toml to find project root
    current = file.resolve()
    while not (current / "pyproject.toml").exists():
        if current.parent == current:
            raise RuntimeError(
                "Could not find project root (pyproject.toml) from {}".format(file)
            )
        current = current.parent
    repo = current
    package = repo / "seqnado"
    # Use test_output as the base for all test artifacts
    test_output = repo / "test_output"
    test_output.mkdir(parents=True, exist_ok=True)
    test_data = test_output / "data"
    test_data.mkdir(parents=True, exist_ok=True)
    workflow = package / "workflow"
    test_profile = workflow / "envs" / "profiles" / "profile_test"
    genome = test_data / "genome"
    genome.mkdir(parents=True, exist_ok=True)
    fastq = test_data / "fastq"
    fastq.mkdir(parents=True, exist_ok=True)
    return TestPaths(
        repo=repo,
        package=package,
        test_dir=test_output,
        test_data=test_data,
        workflow=workflow,
        test_profile=test_profile,
        genome=genome,
        fastq=fastq,
    )


# Unified context class for test session/function-scoped helpers
class TestContext:
    def __init__(self, pytestconfig, tmp_path_factory):
        self.pytestconfig = pytestconfig
        self.tmp_path_factory = tmp_path_factory

        self.test_paths = make_test_paths(Path(__file__).resolve())


    @property
    def cores(self):
        return int(self.pytestconfig.getoption("--cores"))

    @property
    def selected_assays(self):
        return parse_assays(self.pytestconfig)

    def assay_type(self, assay: str) -> str:
        """Extract base assay type from assay name."""
        return re.sub(r"(.*)\-.*", r"\1", assay)

    def run_directory(self, assay: str) -> Path:
        base_temp = self.tmp_path_factory.getbasetemp()
        run_dir = base_temp / assay
        run_dir.mkdir(exist_ok=True)
        return run_dir

    def plot_bed(self, test_data_path: Path):
        return test_data_path / "plotting_coordinates.bed"
