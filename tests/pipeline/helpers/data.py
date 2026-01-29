"""Test data management for SeqNado pipeline tests."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import ClassVar

import requests
from pydantic import BaseModel, field_validator

from .utils import download_with_retry, extract_tar, get_fastq_pattern


class GenomeResources(BaseModel):
    """Validated genome resources for an assay."""

    star_index: Path | None = None
    bt2_index: Path | None = None
    chromsizes: Path | None = None
    gtf: Path | None = None
    blacklist: Path | None = None
    genes_bed: Path | None = None
    fasta: Path | None = None
    viewpoints: Path | None = None

    # Configuration constants
    REFERENCE_URL: ClassVar[str] = (
        "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/seqnado/test_data"
    )

    GENOME_TEMPLATES: ClassVar[dict] = {
        "default": {
            "bt2_index": ("bt2_chr21_dm6_chr2L", "bt2_chr21_dm6_chr2L", "hg38"),
            "chromosome_sizes": "chr21.fa.fai",
            "gtf": "chr21.gtf",
            "blacklist": "hg38-blacklist.v2.bed.gz",
            "genes": "hg38_genes.bed",
        },
        "rna": {
            "star_index": "STAR_chr21_rna_spikein",
            "bt2_index": ("bt2_chr21_dm6_chr2L", "bt2_chr21_dm6_chr2L", "hg38"),
            "chromosome_sizes": "chr21.fa.fai",
            "gtf": "chr21_rna_spikein.gtf",
            "blacklist": "hg38-blacklist.v2.bed.gz",
            "genes": "hg38_genes.bed",
        },
        "meth": {
            "bt2_index": ("bt2_chr21_meth", "chr21_meth", "hg38_meth"),
            "chromosome_sizes": "chr21.fa.fai",
            "gtf": "chr21.gtf",
            "blacklist": "hg38-blacklist.v2.bed.gz",
            "genes": "hg38_genes.bed",
            "fasta": "chr21_meth.fa",
        },
        "crispr": {
            "bt2_index": (
                "bt2_TKOv3_guides_chr21",
                "TKOv3_guides_chr21",
                "TKOv3_guides",
            ),
            "chromosome_sizes": "TKOv3_guides_chr21.fasta.fai",
            "gtf": "TKOv3_guides_chr21.saf",
            "blacklist": "hg38-blacklist.v2.bed.gz",
            "genes": "hg38_genes.bed",
        },
    }

    ASSAY_TEMPLATE_MAP: ClassVar[dict] = {
        "atac": "default",
        "chip": "default",
        "cat": "default",
        "rna": "rna",
        "meth": "meth",
        "snp": "meth",
        "mcc": "default",
        "crispr": "crispr",
    }

    EXTRA_FILES: ClassVar[dict] = {
        "meth": ["chr21_meth.fa.fai"],
        "snp": ["chr21_meth.fa.fai"],
        "mcc": ["chr21_meth.fa", "chr21_meth.fa.fai", "mcc_viewpoints.bed"],
    }

    @field_validator("*", mode="before")
    @classmethod
    def validate_paths(cls, v):
        """Convert string paths to Path objects. Existence is checked at usage time."""
        if v is None:
            return None
        if isinstance(v, (str, Path)):
            return Path(v)
        return v

    def validate_for_assay(self, assay: str) -> None:
        """Validate that required fields for this assay are present."""
        assay_lower = assay.lower()

        # RNA requires star_index and RNA-specific GTF
        if "rna" in assay_lower:
            if self.star_index is None:
                raise FileNotFoundError(
                    f"STAR index required for {assay} but not found"
                )
            if self.gtf is None:
                raise FileNotFoundError(f"GTF required for {assay} but not found")

        # ATAC, ChIP, CAT require bowtie2 index
        elif any(x in assay_lower for x in ["atac", "chip", "cat", "mcc"]):
            if self.bt2_index is None:
                raise FileNotFoundError(
                    f"Bowtie2 index required for {assay} but not found"
                )

        # Methylation and SNP require bowtie2 and fasta
        elif any(x in assay_lower for x in ["meth", "snp"]):
            if self.bt2_index is None:
                raise FileNotFoundError(
                    f"Bowtie2 index required for {assay} but not found"
                )
            if self.fasta is None:
                raise FileNotFoundError(f"FASTA required for {assay} but not found")

        # CRISPR requires bowtie2
        elif "crispr" in assay_lower:
            if self.bt2_index is None:
                raise FileNotFoundError(
                    f"Bowtie2 index required for {assay} but not found"
                )

    def write_config(self, config_file: Path, assay: str | None = None) -> None:
        """Write genome_config.json from this instance.

        Args:
            config_file: Path to write genome_config.json
            assay: Assay key for the config (defaults to 'hg38')
        """
        config_file.parent.mkdir(parents=True, exist_ok=True)

        existing = {}
        if config_file.exists():
            with open(config_file) as f:
                existing = json.load(f)

        genome_key = assay if assay else "hg38"
        existing[genome_key] = {
            "star_index": str(self.star_index) if self.star_index else None,
            "bt2_index": str(self.bt2_index) if self.bt2_index else None,
            "chromosome_sizes": str(self.chromsizes) if self.chromsizes else None,
            "gtf": str(self.gtf) if self.gtf else None,
            "blacklist": str(self.blacklist) if self.blacklist else None,
            "genes": str(self.genes_bed) if self.genes_bed else None,
            "fasta": str(self.fasta) if self.fasta else None,
        }

        with open(config_file, "w") as f:
            json.dump(existing, f, indent=2)

    @classmethod
    def download(cls, genome_path: Path, assay: str) -> GenomeResources:
        """Download genome resources for an assay.

        Args:
            genome_path: Path to store genome files
            assay: Assay type (e.g., 'rna', 'chip', 'atac')

        Returns:
            GenomeResources: Validated genome resources

        Raises:
            FileNotFoundError: If required resources cannot be obtained
        """
        genome_path.mkdir(parents=True, exist_ok=True)

        # Look up template for this assay
        assay_key = next(
            (k for k in cls.ASSAY_TEMPLATE_MAP if k in assay.lower()), "chip"
        )
        template_name = cls.ASSAY_TEMPLATE_MAP.get(assay_key, "default")
        template = cls.GENOME_TEMPLATES[template_name]
        ref_url = f"{cls.REFERENCE_URL}/reference"

        # Use existing factory method
        return cls.from_template(genome_path, template, ref_url, assay)

    @classmethod
    def from_template(
        cls, genome_path: Path, template: dict, ref_url: str, assay: str
    ) -> GenomeResources:
        """Factory method to download resources from template and return validated instance.

        Args:
            genome_path: Path to store genome files
            template: Template dict with required files
            ref_url: Base URL for remote files
            assay: Assay type for validation

        Returns:
            GenomeResources: Validated genome resources

        Raises:
            FileNotFoundError: If required resources cannot be obtained
        """
        config: dict[str, Path | None] = {
            "star_index": None,
            "bt2_index": None,
            "chromsizes": None,
            "gtf": None,
            "blacklist": None,
            "genes_bed": None,
            "fasta": None,
        }

        for key, value in template.items():
            if key == "bt2_index":
                bt2_dir, bt2_prefix, db_name = value
                config["bt2_index"] = cls._ensure_bt2_index(
                    genome_path, bt2_dir, bt2_prefix, db_name, ref_url
                )
            elif key == "star_index":
                config["star_index"] = cls._ensure_star_index(genome_path, ref_url)
            elif key == "chromosome_sizes":
                config["chromsizes"] = cls._ensure_chromsizes(
                    genome_path, value, ref_url
                )
            elif key == "blacklist":
                config["blacklist"] = cls._download(
                    genome_path,
                    value,
                    "https://github.com/Boyle-Lab/Blacklist/raw/master/lists",
                    required=False,  # Blacklist is optional
                )
            elif key == "gtf":
                config["gtf"] = cls._download(
                    genome_path, value, ref_url, required=True
                )
            elif key == "genes":
                config["genes_bed"] = cls._download(
                    genome_path, value, ref_url, required=False
                )
            elif key == "fasta":
                config["fasta"] = cls._download(
                    genome_path, value, ref_url, required=False
                )

        # Download extra files not in config
        assay_key = next(
            (k for k in cls.ASSAY_TEMPLATE_MAP if k in assay.lower()), "chip"
        )
        for filename in cls.EXTRA_FILES.get(assay_key, []):
            cls._download(genome_path, filename, ref_url, required=False)

        # Always download plotting_coordinates
        cls._download(genome_path, "plotting_coordinates.bed", ref_url, required=False)

        # Create and validate instance
        resources = cls(**config)
        resources.validate_for_assay(assay)
        return resources

    @staticmethod
    def _file_exists_remote(url: str) -> bool:
        """Check if a remote file exists without downloading it."""
        try:
            response = requests.head(url, timeout=5, allow_redirects=True)
            return response.status_code < 400
        except Exception:
            return False

    @staticmethod
    def _download(
        genome_path: Path,
        filename: str,
        base_url: str,
        remote_name: str | None = None,
        required: bool = True,
    ) -> Path | None:
        """Download a file if not present. Returns None if optional file is unavailable.

        Args:
            genome_path: Path to genome directory
            filename: Name of file to download
            base_url: Base URL for remote file
            remote_name: Optional different name on remote
            required: If True, raises error if file is unavailable; if False, returns None

        Raises:
            FileNotFoundError: If required file doesn't exist locally or remotely
        """
        dest = genome_path / filename
        if dest.exists():
            return dest

        remote_url = f"{base_url}/{remote_name or filename}"
        if not GenomeResources._file_exists_remote(remote_url):
            msg = f"Remote file not found: {remote_url}"
            if required:
                raise FileNotFoundError(msg)
            print(f"[WARNING] {msg}")
            return None

        try:
            download_with_retry(remote_url, dest)
            return dest
        except Exception as e:
            msg = f"Could not download {filename}: {e}"
            if required:
                raise FileNotFoundError(msg)
            print(f"[WARNING] {msg}")
            return None

    @staticmethod
    def _ensure_star_index(genome_path: Path, ref_url: str) -> Path | None:
        """Download and extract STAR index if not present. Returns None if unavailable."""
        dest = genome_path / "STAR_chr21_rna_spikein"
        if dest.exists():
            return dest

        tar_path = genome_path / "STAR_chr21_rna_spikein.tar.gz"
        try:
            download_with_retry(f"{ref_url}/{tar_path.name}", tar_path)
            extract_tar(tar_path, genome_path)
            tar_path.unlink(missing_ok=True)
        except Exception as e:
            print(f"[WARNING] Could not download STAR index: {e}")

        # Only return if it actually exists
        return dest if dest.exists() else None

    @staticmethod
    def _ensure_bt2_index(
        genome_path: Path, bt2_dir: str, bt2_prefix: str, db_name: str, ref_url: str
    ) -> Path | None:
        """Download and extract Bowtie2 index if not present. Returns None if unavailable."""
        index_dir = genome_path / bt2_dir
        dest = index_dir / bt2_prefix

        if list(index_dir.glob(f"{bt2_prefix}.*.bt2*")):
            return dest

        tar_path = genome_path / f"{bt2_dir}.tar.gz"
        try:
            download_with_retry(f"{ref_url}/{tar_path.name}", tar_path)
            index_dir.mkdir(parents=True, exist_ok=True)
            extract_tar(tar_path, index_dir, flatten=False)
            tar_path.unlink(missing_ok=True)

            # Flatten nested directories if needed
            if not list(index_dir.glob(f"{bt2_prefix}.*.bt2*")):
                for nested in [index_dir / bt2_prefix, index_dir / bt2_dir]:
                    if nested.is_dir():
                        source = (
                            nested / bt2_prefix
                            if (nested / bt2_prefix).is_dir()
                            else nested
                        )
                        for f in source.glob("*.bt2*"):
                            f.rename(index_dir / f.name)
                        try:
                            (nested / bt2_prefix).rmdir()
                        except (OSError, FileNotFoundError):
                            pass
                        try:
                            nested.rmdir()
                        except (OSError, FileNotFoundError):
                            pass
        except Exception as e:
            print(f"[WARNING] Could not download BT2 index: {e}")
            tar_path.unlink(missing_ok=True)

        # Only return if files actually exist
        if list(index_dir.glob(f"{bt2_prefix}.*.bt2*")):
            GenomeResources._update_fastq_screen_config(genome_path, dest, db_name)
            return dest

        return None

    @staticmethod
    def _update_fastq_screen_config(
        genome_path: Path, bt2_index_path: Path, db_name: str
    ) -> None:
        """Append genome entry to fastq_screen.conf.

        Args:
            genome_path: Path to genome directory
            bt2_index_path: Path to the bowtie2 index
            db_name: Database name from template (e.g., 'hg38', 'hg38_meth', 'TKOv3_guides')
        """
        config_path = genome_path / "fastq_screen.conf"

        # Determine database name based on the index being used
        if "meth" in db_name.lower():
            fs_db_name = "hg38_meth"
        elif "TKOv3" in db_name:
            fs_db_name = "TKOv3_guides"
        else:
            fs_db_name = "hg38"

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
        if fs_db_name not in existing_entries:
            with open(config_path, "a") as f:
                f.write(f"DATABASE\t{fs_db_name}\t{bt2_index_path}\n")

    @staticmethod
    def _ensure_chromsizes(
        genome_path: Path, chromsizes_file: str, ref_url: str
    ) -> Path:
        """Download chromosome sizes if not present."""
        dest = genome_path / chromsizes_file

        if chromsizes_file != "chr21.fa.fai":
            return GenomeResources._download(genome_path, chromsizes_file, ref_url)

        if (genome_path / "chr21_rename.fa.fai").exists():
            return dest

        download_with_retry(f"{ref_url}/{dest.name}", dest)

        # Append dm6 chromosome sizes
        for attempt in range(3):
            try:
                r = requests.get(
                    "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes",
                    stream=True,
                    timeout=30,
                )
                r.raise_for_status()
                with open(dest, "ab") as f:
                    for line in r.iter_lines():
                        f.write(b"dm6_" + line + b"\n")
                return dest
            except requests.RequestException:
                if attempt < 2:
                    time.sleep(2**attempt)
                else:
                    raise
        return dest


class FastqFiles(BaseModel):
    """Manages FASTQ file downloads and organization by assay."""

    fastq_dir: Path
    files: dict[str, list[Path]] = {}

    @classmethod
    def download(cls, fastq_dir: Path, selected_assays: list[str]) -> FastqFiles:
        """Download FASTQ files and organize them by assay.

        Args:
            fastq_dir: Directory to store FASTQ files
            selected_assays: List of assay types to get FASTQs for

        Returns:
            FastqFiles: Instance with files organized by assay

        Raises:
            RuntimeError: If FASTQ extraction fails
        """
        fastq_dir.mkdir(parents=True, exist_ok=True)

        if not list(fastq_dir.glob("*.fastq.gz")):
            tar_path = fastq_dir.parent / "fastq.tar.gz"
            download_with_retry(
                f"{GenomeResources.REFERENCE_URL}/fastq.tar.gz", tar_path
            )
            extract_tar(tar_path, fastq_dir.parent, flatten=False)

            nested = fastq_dir / "fastq"
            if nested.is_dir():
                for f in nested.glob("*.fastq.gz"):
                    f.rename(fastq_dir / f.name)
                nested.rmdir()

            tar_path.unlink(missing_ok=True)

            if not list(fastq_dir.glob("*.fastq.gz")):
                raise RuntimeError(f"FASTQ extraction failed - no files in {fastq_dir}")

        files = {
            a: sorted(fastq_dir.glob(get_fastq_pattern(a))) for a in selected_assays
        }
        return cls(fastq_dir=fastq_dir, files=files)
