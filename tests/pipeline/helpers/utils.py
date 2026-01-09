import json
import tarfile
import time
from pathlib import Path

import requests


def download_with_retry(url: str, dest: Path, max_retries: int = 3, timeout: int = 30):
    """Download a file with retry logic."""
    for attempt in range(max_retries):
        try:
            dest.parent.mkdir(parents=True, exist_ok=True)
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
    """Extract a tar archive, optionally flattening nested directories."""
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


def get_fastq_pattern(assay: str) -> str:
    """Get the FASTQ file pattern for a given assay."""
    patterns = {
        "atac": "atac_*.fastq.gz",
        "cat": "chip-rx_*.fastq.gz",
        "chip-rx": "chip-rx_*.fastq.gz",
        "chip": "chip-rx_*.fastq.gz",
        "crispr": "crispr_*.fastq.gz",
        "mcc": "mcc_*.fastq.gz",
        "meth": "meth-*.fastq.gz",
        "rna-rx": "rna-spikein-*.fastq.gz",
        "rna": "rna_*.fastq.gz",
        "snp": "snp_*.fastq.gz",
    }
    pattern = patterns.get(assay)
    if not pattern:
        raise ValueError(f"Unsupported assay: {assay}")
    return pattern


def setup_genome_config(
    genome_config_file: Path,
    star_index: Path,
    bt2_index: Path,
    chromsizes: Path,
    gtf: Path,
    blacklist: Path,
    genes_bed: Path,
    fasta: Path,
    assay: str = None,
):
    """Write genome_config.json for seqnado tests.

    If assay is provided, creates an assay-specific genome entry (e.g., "meth").
    Otherwise, creates the default "hg38" entry.

    For Multiomic tests, this function can be called multiple times to add
    assay-specific configurations to the same file.
    """
    genome_config_file.parent.mkdir(parents=True, exist_ok=True)

    # Load existing config if it exists, but filter to only keep test-specific entries
    # Valid test genome keys are: assay names (atac, chip, meth, rna, snp, etc.) or "hg38"
    valid_test_genomes = ["hg38", "atac", "chip", "cat", "chip-rx", "meth", "rna", "rna-rx", "snp", "mcc"]

    if genome_config_file.exists():
        with open(genome_config_file, "r") as f:
            all_config_data = json.load(f)
        # Filter to keep only test-specific genome configs
        config_data = {k: v for k, v in all_config_data.items() if k in valid_test_genomes}
    else:
        config_data = {}

    # Determine the genome key - use assay name if provided, otherwise "hg38"
    genome_key = assay if assay else "hg38"

    # Add or update the genome configuration for this assay
    config_data[genome_key] = {
        "star_index": str(star_index) if star_index else None,
        "bt2_index": str(bt2_index) if bt2_index else None,
        "chromosome_sizes": str(chromsizes) if chromsizes else None,
        "gtf": str(gtf) if gtf else None,
        "blacklist": str(blacklist) if blacklist else None,
        "genes": str(genes_bed) if genes_bed else None,
        "fasta": str(fasta) if fasta else None,
    }

    with open(genome_config_file, "w") as f:
        json.dump(config_data, f, indent=2)
