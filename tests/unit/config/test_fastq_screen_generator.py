"""Tests for fastq_screen_generator module."""

import json
from pathlib import Path

import pytest

from seqnado.config.configs import BowtieIndex, GenomeConfig
from seqnado.config.fastq_screen_generator import (
    generate_fastq_screen_config,
    load_genome_configs_for_fastqscreen,
)


@pytest.fixture
def test_genome_dir(tmp_path):
    """Create temporary test genome data directory."""
    genome_dir = tmp_path / "genome"
    genome_dir.mkdir(parents=True)
    
    # Create Bowtie2 index directory
    bt2_dir = genome_dir / "bt2_chr21_dm6_chr2L"
    bt2_dir.mkdir(parents=True)
    for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
        (bt2_dir / f"bt2_chr21_dm6_chr2L{suffix}").touch()
    
    # Create fasta index file
    (genome_dir / "chr21.fa.fai").write_text("chr21\t48129895\t4\t50\t51\n")
    
    return genome_dir


@pytest.fixture
def bt2_index_path(test_genome_dir):
    """Get the Bowtie2 index path."""
    return test_genome_dir / "bt2_chr21_dm6_chr2L" / "bt2_chr21_dm6_chr2L"


@pytest.fixture
def mock_genome_configs(tmp_path):
    """Create mock genome configs with temporary file structure."""
    # Create the directory structure that BowtieIndex expects
    bt2_dir = tmp_path / "bt2_chr21_dm6_chr2L"
    bt2_dir.mkdir(parents=True, exist_ok=True)
    
    # Touch the index files that would exist
    for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
        (bt2_dir / f"bt2_chr21_dm6_chr2L{suffix}").touch()
    
    # Create fasta index file
    fasta_file = tmp_path / "chr21.fa.fai"
    fasta_file.touch()
    
    bt2_prefix = bt2_dir / "bt2_chr21_dm6_chr2L"
    
    return {
        "hg38": GenomeConfig(
            name="hg38",
            fasta=str(fasta_file),
            index=BowtieIndex(prefix=str(bt2_prefix)),
        ),
        "mm39": GenomeConfig(
            name="mm39",
            fasta=str(fasta_file),
            index=BowtieIndex(prefix=str(bt2_prefix)),
        ),
        "mm10": GenomeConfig(  # Older version, should be skipped
            name="mm10",
            fasta=str(fasta_file),
            index=BowtieIndex(prefix=str(bt2_prefix)),
        ),
    }


@pytest.fixture
def temp_output_path(tmp_path):
    """Create a temporary output path."""
    return tmp_path / "fastq_screen.conf"


def test_generate_fastq_screen_config_basic(
    mock_genome_configs, temp_output_path, tmp_path
):
    """Test basic fastq_screen config generation without contaminants."""
    # Generate config without contaminants
    generate_fastq_screen_config(
        genome_configs=mock_genome_configs,
        output_path=temp_output_path,
        threads=4,
        include_contaminants=False,
    )

    # Check file was created
    assert temp_output_path.exists()

    # Read and verify content
    content = temp_output_path.read_text()

    # Check header elements
    assert "FastqScreen configuration file" in content
    assert "THREADS\t\t4" in content

    # Should have Human (hg38) and Mouse (mm39, not mm10)
    assert "Human" in content
    assert "Mouse" in content


def test_generate_fastq_screen_config_with_contaminants(
    mock_genome_configs, temp_output_path, tmp_path
):
    """Test fastq_screen config generation with contaminant databases."""
    contaminant_path = tmp_path / "contaminants"
    contaminant_path.mkdir()

    # Create mock contaminant directories and index files
    ecoli_dir = contaminant_path / "E_coli"
    ecoli_dir.mkdir()
    (ecoli_dir / "Ecoli.1.bt2").touch()

    generate_fastq_screen_config(
        genome_configs=mock_genome_configs,
        output_path=temp_output_path,
        threads=8,
        include_contaminants=True,
        contaminant_base_path=contaminant_path,
    )

    assert temp_output_path.exists()
    content = temp_output_path.read_text()

    # Should have Ecoli contaminant
    assert "Ecoli" in content


def test_generate_fastq_screen_config_custom_threads(
    mock_genome_configs, temp_output_path
):
    """Test custom thread count in config."""
    generate_fastq_screen_config(
        genome_configs=mock_genome_configs,
        output_path=temp_output_path,
        threads=16,
        include_contaminants=False,
    )

    content = temp_output_path.read_text()
    assert "THREADS\t\t16" in content


def test_generate_fastq_screen_config_creates_parent_dir(mock_genome_configs, tmp_path):
    """Test that parent directories are created if they don't exist."""
    nested_path = tmp_path / "nested" / "dir" / "fastq_screen.conf"

    generate_fastq_screen_config(
        genome_configs=mock_genome_configs,
        output_path=nested_path,
        include_contaminants=False,
    )

    assert nested_path.exists()
    assert nested_path.parent.exists()


def test_load_genome_configs_for_fastqscreen_success(tmp_path, monkeypatch, test_genome_dir, bt2_index_path):
    """Test loading genome configs from config file."""
    # Create mock config file
    config_dir = tmp_path / ".config" / "seqnado"
    config_dir.mkdir(parents=True)
    config_file = config_dir / "genome_config.json"

    config_data = {
        "hg38": {
            "fasta": str(test_genome_dir / "chr21.fa.fai"),
            "bt2_index": str(bt2_index_path),
        },
        "mm39": {
            "fasta": str(test_genome_dir / "chr21.fa.fai"),
            "bt2_index": str(bt2_index_path),
        },
    }

    config_file.write_text(json.dumps(config_data))

    # Mock environment variable
    monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))

    # Load configs
    configs = load_genome_configs_for_fastqscreen()

    assert "hg38" in configs
    assert "mm39" in configs
    assert configs["hg38"].name == "hg38"
    assert str(bt2_index_path) in configs["hg38"].index.prefix


def test_load_genome_configs_for_fastqscreen_missing_file(tmp_path, monkeypatch):
    """Test error when config file doesn't exist."""
    monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))

    with pytest.raises(FileNotFoundError, match="Genome config not found"):
        load_genome_configs_for_fastqscreen()


def test_generate_fastq_screen_config_skips_combined_genomes(tmp_path, test_genome_dir, bt2_index_path):
    """Test that combined genomes (e.g., dm6_mm39) are skipped."""
    configs = {
        "dm6": GenomeConfig(
            name="dm6",
            fasta=str(test_genome_dir / "chr21.fa.fai"),
            index=BowtieIndex(prefix=str(bt2_index_path)),
        ),
        "dm6_mm39": GenomeConfig(  # Combined genome - should be skipped
            name="dm6_mm39",
            fasta=str(test_genome_dir / "chr21.fa.fai"),
            index=BowtieIndex(prefix=str(bt2_index_path)),
        ),
    }

    output_path = tmp_path / "test_fastq_screen.conf"

    generate_fastq_screen_config(
        genome_configs=configs,
        output_path=output_path,
        include_contaminants=False,
    )

    content = output_path.read_text()

    # Should have dm6
    assert "Drosophila" in content or "dm6" in content
    # Combined genome name should not appear as database
    assert "dm6_mm39" not in content


def test_generate_fastq_screen_config_missing_index_files(tmp_path, capsys, test_genome_dir):
    """Test warning when index files don't exist."""
    # Create directory with index files first
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    # Create minimal index files to pass validation
    (empty_dir / "test.1.bt2").touch()

    # Create the config
    config = GenomeConfig(
        name="hg38",
        fasta=str(test_genome_dir / "chr21.fa.fai"),
        index=BowtieIndex(prefix=str(empty_dir / "test")),
    )

    # Now remove the index files to trigger the warning
    for f in empty_dir.glob("*.bt2"):
        f.unlink()

    configs = {"hg38": config}

    output_path = tmp_path / "test.conf"
    generate_fastq_screen_config(
        genome_configs=configs,
        output_path=output_path,
        include_contaminants=False,
    )

    # Check warning was printed
    captured = capsys.readouterr()
    assert "Warning: No bowtie2 index files found" in captured.out


def test_generate_fastq_screen_config_missing_index_directory(tmp_path, capsys, test_genome_dir):
    """Test warning when index directory doesn't exist."""
    # Create directory and index files first
    temp_dir = tmp_path / "temp"
    temp_dir.mkdir()
    (temp_dir / "hg38.1.bt2").touch()

    config = GenomeConfig(
        name="hg38",
        fasta=str(test_genome_dir / "chr21.fa.fai"),
        index=BowtieIndex(prefix=str(temp_dir / "hg38")),
    )

    # Remove the directory to trigger the warning
    for f in temp_dir.glob("*"):
        f.unlink()
    temp_dir.rmdir()

    configs = {"hg38": config}

    output_path = tmp_path / "test.conf"
    generate_fastq_screen_config(
        genome_configs=configs,
        output_path=output_path,
        include_contaminants=False,
    )

    # Check warning was printed
    captured = capsys.readouterr()
    assert "Warning: Index path does not exist" in captured.out


def test_generate_fastq_screen_config_unknown_organism(tmp_path, test_genome_dir, bt2_index_path):
    """Test display name for unknown organism (custom genome)."""
    configs = {
        "custom_genome": GenomeConfig(
            name="custom_genome",
            fasta=str(test_genome_dir / "chr21.fa.fai"),
            index=BowtieIndex(prefix=str(bt2_index_path)),
        ),
    }

    output_path = tmp_path / "test.conf"
    generate_fastq_screen_config(
        genome_configs=configs,
        output_path=output_path,
        include_contaminants=False,
    )

    content = output_path.read_text()
    # Should convert underscores and title case
    assert "Custom Genome" in content


def test_generate_fastq_screen_config_contaminants_no_path(tmp_path, capsys, test_genome_dir, bt2_index_path):
    """Test warning when include_contaminants=True but no path provided."""
    configs = {
        "hg38": GenomeConfig(
            name="hg38",
            fasta=str(test_genome_dir / "chr21.fa.fai"),
            index=BowtieIndex(prefix=str(bt2_index_path)),
        ),
    }

    output_path = tmp_path / "test.conf"
    generate_fastq_screen_config(
        genome_configs=configs,
        output_path=output_path,
        include_contaminants=True,
        contaminant_base_path=None,
    )

    captured = capsys.readouterr()
    assert (
        "include_contaminants=True but no contaminant_base_path provided"
        in captured.out
    )
    assert "Skipping contaminant databases" in captured.out


def test_generate_fastq_screen_config_contaminants_path_not_exist(tmp_path, capsys, test_genome_dir, bt2_index_path):
    """Test warning when contaminant path doesn't exist."""
    nonexistent_path = tmp_path / "nonexistent_contaminants"
    
    configs = {
        "hg38": GenomeConfig(
            name="hg38",
            fasta=str(test_genome_dir / "chr21.fa.fai"),
            index=BowtieIndex(prefix=str(bt2_index_path)),
        ),
    }
    
    output_path = tmp_path / "test.conf"
    generate_fastq_screen_config(
        genome_configs=configs,
        output_path=output_path,
        include_contaminants=True,
        contaminant_base_path=nonexistent_path,
    )
    
    captured = capsys.readouterr()
    assert "Contaminant path does not exist" in captured.out
    assert "Skipping contaminant databases" in captured.out


def test_generate_fastq_screen_config_version_selection(test_genome_dir, bt2_index_path):
    """Test that newer genome version is selected over older."""
    configs = {
        "mm10": GenomeConfig(
            name="mm10",
            fasta=str(test_genome_dir / "chr21.fa.fai"),
            index=BowtieIndex(prefix=str(bt2_index_path)),
        ),
        "mm39": GenomeConfig(
            name="mm39",
            fasta=str(test_genome_dir / "chr21.fa.fai"),
            index=BowtieIndex(prefix=str(bt2_index_path)),
        ),
        "hg19": GenomeConfig(
            name="hg19",
            fasta=str(test_genome_dir / "chr21.fa.fai"),
            index=BowtieIndex(prefix=str(bt2_index_path)),
        ),
        "hg38": GenomeConfig(
            name="hg38",
            fasta=str(test_genome_dir / "chr21.fa.fai"),
            index=BowtieIndex(prefix=str(bt2_index_path)),
        ),
    }
    
    output_path = Path("/tmp/test_version.conf")
    generate_fastq_screen_config(
        genome_configs=configs,
        output_path=output_path,
        include_contaminants=False,
    )
    
    content = output_path.read_text()
    
    # Should have hg38 and mm39 (newer versions)
    assert "hg38" in content or "Human" in content
    assert "mm39" in content or "Mouse" in content
    # Should NOT have hg19 or mm10 (they're replaced by newer versions)
    assert "hg19" not in content
    assert "mm10" not in content
    
    output_path.unlink()


def test_generate_fastq_screen_config_dict_index(bt2_index_path):
    """Test handling of dict-based index (alternative to object)."""
    # Create a GenomeConfig with dict-like index access
    from unittest.mock import Mock
    
    mock_config = Mock(spec=GenomeConfig)
    mock_config.name = "test_genome"
    # Simulate dict index (no 'prefix' attribute)
    mock_config.index = {"prefix": str(bt2_index_path)}
    
    configs = {"test_genome": mock_config}
    
    output_path = Path("/tmp/test_dict_index.conf")
    generate_fastq_screen_config(
        genome_configs=configs,
        output_path=output_path,
        include_contaminants=False,
    )
    
    content = output_path.read_text()
    # Should have the genome with title-cased name
    assert "Test Genome" in content
    
    output_path.unlink()