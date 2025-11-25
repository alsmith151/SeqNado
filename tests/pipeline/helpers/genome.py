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
        - fasta: FASTA file (if methylation assay)
        - fasta_fai: FASTA index (if methylation assay)
        - viewpoints: MCC viewpoints file (if MCC assay)
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

    # Genome FASTA files (for methylation)
    if "meth" in assay.lower():
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

    if dest.exists():
        return dest

    url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/{suffix}"
    tar_path = genome_path / suffix

    download_with_retry(url, tar_path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    extract_tar(tar_path, dest.parent)

    # Handle nested directory structure
    nested = dest.parent / dest.name
    if nested.exists() and nested.is_dir():
        for f in nested.glob("*.bt2*"):
            f.rename(dest.parent / f.name)
        try:
            nested.rmdir()
        except OSError:
            pass

    os.remove(tar_path)
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


def fill_fastq_screen_config(config_path, genome_bt2_index):
    """
    Fill the fastq_screen.conf file with the provided genome Bowtie2 index.

    Args:
        config_path (str): Path to the fastq_screen.conf file.
        genome_bt2_index (str): Path to the Bowtie2 index for the genome.
    """
    with open(config_path, "a") as f:
        f.write("# Genome configuration\n")
        f.write(f"DATABASE\thg38\t{genome_bt2_index}\n")
