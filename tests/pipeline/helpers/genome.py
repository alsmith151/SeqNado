"""Genome resource download helpers for pipeline tests."""
import os
import time
from pathlib import Path

import requests

from .utils import download_with_retry, extract_tar


def ensure_genome_resources(genome_path: Path, assay: str) -> dict[str, Path]:
    """
    Ensure all required genome resources are downloaded for the given assay.

    Args:
        genome_path: Base directory for genome data
        assay: Assay type (e.g., 'chip', 'rna', 'atac', 'meth', 'mcc')

    Returns:
        Dictionary with paths to all genome resources:
        - star_index: STAR index (if RNA assay)
        - bt2_index: Bowtie2 index
        - chromsizes: Chromosome sizes file
        - gtf: GTF annotation file
        - blacklist: Blacklist BED file
        - fasta: FASTA file (if methylation assay or snp)
        - fasta_fai: FASTA index (if methylation assay or snp)
        - viewpoints: MCC viewpoints file (if MCC assay)
        -
    """
    genome_path.mkdir(parents=True, exist_ok=True)
    resources = {}

    # STAR index for RNA assays
    if "rna" in assay.lower():
        resources["star_index"] = _ensure_star_index(genome_path)

    # Bowtie2 index
    resources["bt2_index"] = _ensure_bt2_index(genome_path, assay)

    # Chromsizes
    resources["chromsizes"] = _ensure_chromsizes(genome_path)

    # GTF
    resources["gtf"] = _ensure_gtf(genome_path, assay)

    # Blacklist
    resources["blacklist"] = _ensure_blacklist(genome_path)

    # Genome FASTA files (for methylation or SNP assays)
    if any(x in assay.lower() for x in ["meth", "snp",'mcc']):
        fasta, fasta_fai = _ensure_genome_fasta(genome_path)
        resources["fasta"] = fasta
        resources["fasta_fai"] = fasta_fai

    # MCC viewpoints
    if assay.lower() == "mcc":
        resources["viewpoints"] = _ensure_mcc_viewpoints(genome_path)

    return resources


def _ensure_star_index(genome_path: Path) -> Path:
    """Download and extract STAR index if not present."""
    dest = genome_path / "STAR_chr21_rna_spikein"
    if dest.exists():
        return dest

    suffix = "STAR_chr21_rna_spikein.tar.gz"
    url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/{suffix}"
    tar_path = genome_path / suffix

    download_with_retry(url, tar_path)
    extract_tar(tar_path, genome_path)

    # Handle nested directory structure
    nested = dest / dest.name
    if nested.exists() and nested.is_dir():
        for f in nested.iterdir():
            f.rename(dest / f.name)
        nested.rmdir()

    os.remove(tar_path)
    return dest


def _ensure_bt2_index(genome_path: Path, assay: str) -> Path:
    """Download and extract Bowtie2 index if not present."""
    if "meth" in assay.lower():
        dest = genome_path / "bt2_chr21_meth" / "chr21_meth"
        suffix = "bt2_chr21_meth.tar.gz"
    else:
        dest = genome_path / "bt2_chr21_dm6_chr2L" / "bt2_chr21_dm6_chr2L"
        suffix = "bt2_chr21_dm6_chr2L.tar.gz"

    # Check if index files already exist (not just the directory)
    index_files = list(dest.parent.glob(f"{dest.name}.*.bt2*"))
    if index_files:
        return dest

    url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/{suffix}"
    tar_path = genome_path / suffix

    download_with_retry(url, tar_path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Extract without flattening, we'll handle the directory structure manually
    extract_tar(tar_path, dest.parent, flatten=False)

    # Handle nested directory structure - tar might contain bt2_chr21_dm6_chr2L/bt2_chr21_dm6_chr2L/*.bt2
    # or bt2_chr21_dm6_chr2L/*.bt2, so we need to check both levels
    nested = dest.parent / dest.name
    if nested.exists() and nested.is_dir():
        # Check if files are in a doubly-nested directory
        doubly_nested = nested / dest.name
        if doubly_nested.exists() and doubly_nested.is_dir():
            # Move files from bt2_chr21_dm6_chr2L/bt2_chr21_dm6_chr2L/ to bt2_chr21_dm6_chr2L/
            for f in doubly_nested.glob("*.bt2*"):
                f.rename(nested / f.name)
            try:
                doubly_nested.rmdir()
            except OSError:
                pass

    os.remove(tar_path)

    # Verify that index files exist
    index_files = list(dest.parent.glob(f"{dest.name}.*.bt2*"))
    if not index_files:
        raise RuntimeError(
            f"Bowtie2 index extraction failed - no .bt2* files found in {dest.parent}. "
            f"Expected files like {dest.name}.1.bt2, etc."
        )

    # Append genome entry to fastq_screen.conf
    genome_bt2_index = str(dest)
    config_path = genome_path / "fastq_screen.conf"

    # Determine database name based on the index being used
    # Use the index name to create a unique database identifier
    if "meth" in assay.lower():
        db_name = "hg38_meth"
    else:
        db_name = "hg38"

    # Check if this database entry already exists to avoid duplicates
    existing_entries = set()
    if config_path.exists():
        with open(config_path, "r") as f:
            for line in f:
                if line.startswith("DATABASE"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        existing_entries.add(parts[1])

    # Only append if this database name doesn't already exist
    if db_name not in existing_entries:
        with open(config_path, "a") as f:
            f.write(f"DATABASE\t{db_name}\t{genome_bt2_index}\n")

    return dest


def _ensure_chromsizes(genome_path: Path) -> Path:
    """Download chromosome sizes if not present."""
    dest = genome_path / "chr21.fa.fai"

    # Check if already downloaded and processed
    if (genome_path / "chr21_rename.fa.fai").exists():
        return dest

    # Download human chromosome 21 sizes
    url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/chr21.fa.fai"
    download_with_retry(url, dest)

    # Append Drosophila chromosome sizes with dm6_ prefix
    with open(dest, "ab") as f:
        url2 = "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes"
        for attempt in range(3):
            try:
                r = requests.get(url2, stream=True, timeout=30)
                r.raise_for_status()
                for line in r.iter_lines():
                    f.write(b"dm6_" + line + b"\n")
                break
            except requests.RequestException:
                if attempt < 2:
                    time.sleep(2**attempt)
                else:
                    raise

    return dest


def _ensure_gtf(genome_path: Path, assay: str) -> Path:
    """Download GTF annotation file if not present."""
    if "rna" in assay.lower():
        gtf_path = genome_path / "chr21_rna_spikein.gtf"
    else:
        gtf_path = genome_path / "chr21.gtf"

    if gtf_path.exists():
        return gtf_path

    url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/{gtf_path.name}"
    download_with_retry(url, gtf_path)
    return gtf_path


def _ensure_blacklist(genome_path: Path) -> Path:
    """Download ENCODE blacklist file if not present."""
    dest = genome_path / "hg38-blacklist.v2.bed.gz"

    if dest.exists():
        return dest

    url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
    download_with_retry(url, dest)
    return dest


def _ensure_genome_fasta(genome_path: Path) -> tuple[Path, Path]:
    """Download genome FASTA and index files if not present."""
    fasta = genome_path / "chr21_meth.fa"
    fasta_fai = genome_path / "chr21_meth.fa.fai"

    if not fasta.exists():
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/chr21_meth.fa"
        download_with_retry(url, fasta)

    if not fasta_fai.exists():
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/chr21_meth.fa.fai"
        download_with_retry(url, fasta_fai)

    return fasta, fasta_fai


def _ensure_mcc_viewpoints(genome_path: Path) -> Path:
    """Download MCC viewpoints file if not present."""
    dest = genome_path / "mcc_viewpoints.bed"

    if dest.exists():
        return dest

    url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/seqnado_data/test_viewpoints.bed"
    download_with_retry(url, dest)
    return dest
