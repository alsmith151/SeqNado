import os
from pathlib import Path

from .utils import download_with_retry, extract_tar, get_fastq_pattern


def ensure_fastqs_present(
    fastq_dir: Path, selected_assays: list[str]
) -> dict[str, list[Path]]:
    """
    Download FASTQ files and return them indexed by assay.

    Args:
        fastq_dir: Directory where FASTQ files should be stored
        selected_assays: List of assay names to get FASTQs for

    Returns:
        Dict mapping assay name to list of FASTQ file paths
    """
    # Download FASTQs if not present
    fastq_dir.mkdir(parents=True, exist_ok=True)
    existing_fastqs = list(fastq_dir.glob("*.fastq.gz"))

    if not existing_fastqs:
        import tarfile

        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/seqnado_data/fastq.tar.gz"
        tar_path = fastq_dir.parent / "fastq.tar.gz"
        download_with_retry(url, tar_path)

        # Extract tar directly to fastq_dir
        with tarfile.open(tar_path) as tar:
            tar.extractall(path=fastq_dir)

        # Verify extraction worked
        extracted_fastqs = list(fastq_dir.glob("*.fastq.gz"))
        if not extracted_fastqs:
            raise RuntimeError(f"FASTQ extraction failed - no files found in {fastq_dir} after extraction")

        os.remove(tar_path)

    def pick_fastqs(assay_name: str) -> list[Path]:
        pattern = get_fastq_pattern(assay_name)
        
        return sorted(fastq_dir.glob(pattern))

    return {a: pick_fastqs(a) for a in selected_assays}


def ensure_all_genome_data(genome_path: Path, genome_files: dict):
    """
    Download all genome data files and extract indices as needed.

    Args:
        genome_path: Base directory for genome data
        genome_files: Dict mapping filename to download URL
    """
    for fname, url in genome_files.items():
        dest = genome_path / fname
        if not dest.exists():
            download_with_retry(url, dest)

    # Extract STAR index
    star_dir = genome_path / "STAR_chr21_rna_spikein"
    star_tar = genome_path / "STAR_chr21_rna_spikein.tar.gz"
    if star_tar.exists() and not star_dir.exists():
        extract_tar(star_tar, genome_path, flatten=True)

    # Extract Bowtie2 index
    bt2_dir = genome_path / "bt2_chr21_dm6_chr2L" / "bt2_chr21_dm6_chr2L"
    bt2_tar = genome_path / "bt2_chr21_dm6_chr2L.tar.gz"
    if bt2_tar.exists() and not bt2_dir.exists():
        dest_parent = bt2_dir.parent
        dest_parent.mkdir(parents=True, exist_ok=True)
        extract_tar(bt2_tar, dest_parent)
        # Move .bt2* files up if nested
        nested = dest_parent / bt2_dir.name
        if nested.exists() and nested.is_dir():
            for f in nested.glob("*.bt2*"):
                f.rename(dest_parent / f.name)
            try:
                nested.rmdir()
            except Exception as e:
                print(f"[WARNING] Could not remove nested directory {nested}: {e}")

    # Extract Methylation Bowtie2 index
    meth_index_dir = genome_path / "bt2_chr21_meth"
    meth_index_tar = genome_path / "bt2_chr21_meth.tar.gz"
    if meth_index_tar.exists() and not meth_index_dir.exists():
        dest_parent = meth_index_dir.parent
        dest_parent.mkdir(parents=True, exist_ok=True)
        extract_tar(meth_index_tar, dest_parent)
        # Move .bt2* files up if nested
        nested = dest_parent / meth_index_dir.name
        if nested.exists() and nested.is_dir():
            for f in nested.glob("*.bt2*"):
                f.rename(dest_parent / f.name)
            try:
                nested.rmdir()
            except Exception as e:
                print(f"[WARNING] Could not remove nested directory {nested}: {e}")
