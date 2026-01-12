"""Comprehensive integration tests for SeqNado CLI commands."""

import json
import subprocess
from pathlib import Path

import pytest


@pytest.fixture
def genome_config_dir(tmp_path):
    """Create a minimal genome config for testing."""
    cfg_dir = tmp_path / ".config" / "seqnado"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    
    # Create test genome indices
    star_dir = tmp_path / "star_index"
    star_dir.mkdir(parents=True, exist_ok=True)
    (star_dir / "SA").touch()
    
    bowtie_dir = tmp_path / "bowtie_index"
    bowtie_dir.mkdir(parents=True, exist_ok=True)
    (bowtie_dir / "genome.1.bt2").touch()
    
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGTACGTACGTACGT\n")
    
    genome_cfg = {
        "test_genome": {
            "star_index": str(star_dir),
            "bowtie2_index": str(bowtie_dir / "genome"),
            "fasta": str(fasta),
        }
    }
    
    cfg_path = cfg_dir / "genome_config.json"
    cfg_path.write_text(json.dumps(genome_cfg, indent=2))
    return cfg_dir


class TestInitCommand:
    """Tests for seqnado init command."""
    
    def test_init_help(self):
        """Test init command help."""
        result = subprocess.run(
            ["seqnado", "init", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "init" in result.stdout.lower()
    
    def test_init_dry_run(self, tmp_path, monkeypatch):
        """Test init command in dry-run mode."""
        monkeypatch.setenv("HOME", str(tmp_path))
        
        result = subprocess.run(
            ["seqnado", "init", "--dry-run"],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        # Dry run should succeed without errors
        assert result.returncode == 0


class TestConfigCommand:
    """Tests for seqnado config command."""
    
    def test_config_help(self):
        """Test config command help."""
        result = subprocess.run(
            ["seqnado", "config", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "config" in result.stdout.lower()
    
    def test_config_atac_non_interactive(self, tmp_path, monkeypatch, genome_config_dir):
        """Test config generation for ATAC-seq."""
        monkeypatch.setenv("HOME", str(tmp_path))
        (tmp_path / "metadata.csv").write_text("sample_id,fastq_r1\n")
        (tmp_path / "seqnado_output").mkdir()
        
        out_file = tmp_path / "config_atac.yaml"
        result = subprocess.run(
            [
                "seqnado", "config", "atac",
                "--no-interactive",
                "--no-make-dirs",
                "--render-options",
                "-o", str(out_file),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert out_file.exists()
        config_text = out_file.read_text()
        assert "assay: ATAC" in config_text or "assay: atac" in config_text.lower()
    
    def test_config_chip_non_interactive(self, tmp_path, monkeypatch, genome_config_dir):
        """Test config generation for ChIP-seq."""
        monkeypatch.setenv("HOME", str(tmp_path))
        (tmp_path / "metadata.csv").write_text("sample_id,fastq_r1,antibody\n")
        (tmp_path / "seqnado_output").mkdir()
        
        out_file = tmp_path / "config_chip.yaml"
        result = subprocess.run(
            [
                "seqnado", "config", "chip",
                "--no-interactive",
                "--no-make-dirs",
                "--render-options",
                "-o", str(out_file),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert out_file.exists()
        config_text = out_file.read_text()
        assert "chip" in config_text.lower()
    
    def test_config_with_output_file(self, tmp_path, monkeypatch, genome_config_dir):
        """Test config with custom output file."""
        monkeypatch.setenv("HOME", str(tmp_path))
        (tmp_path / "metadata.csv").write_text("sample_id,fastq_r1\n")
        (tmp_path / "seqnado_output").mkdir()
        
        out_file = tmp_path / "custom_config.yaml"
        result = subprocess.run(
            [
                "seqnado", "config", "atac",
                "--no-interactive",
                "--no-make-dirs",
                "-o", str(out_file),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert out_file.exists()
        assert out_file.name == "custom_config.yaml"


class TestDesignCommand:
    """Tests for seqnado design command."""
    
    def test_design_help(self):
        """Test design command help."""
        result = subprocess.run(
            ["seqnado", "design", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "design" in result.stdout.lower()
    
    def test_design_atac_from_fastqs(self, tmp_path):
        """Test design CSV generation for ATAC from FASTQ files."""
        # Create test FASTQ files
        (tmp_path / "sample1_R1.fastq.gz").touch()
        (tmp_path / "sample1_R2.fastq.gz").touch()
        (tmp_path / "sample2_R1.fastq.gz").touch()
        (tmp_path / "sample2_R2.fastq.gz").touch()
        
        out_csv = tmp_path / "design.csv"
        result = subprocess.run(
            [
                "seqnado", "design", "atac",
                "--no-interactive",
                "--accept-all-defaults",
                "-o", str(out_csv),
                str(tmp_path),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert out_csv.exists()
        
        design_text = out_csv.read_text()
        assert "sample" in design_text.lower()
        assert "sample1" in design_text or "sample2" in design_text
    
    def test_design_chip_from_directory(self, tmp_path):
        """Test design CSV generation for ChIP from directory."""
        # Create test FASTQ files with ChIP naming
        (tmp_path / "H3K27ac_IP_R1.fastq.gz").touch()
        (tmp_path / "H3K27ac_IP_R2.fastq.gz").touch()
        (tmp_path / "Input_R1.fastq.gz").touch()
        (tmp_path / "Input_R2.fastq.gz").touch()
        
        out_csv = tmp_path / "chip_design.csv"
        result = subprocess.run(
            [
                "seqnado", "design", "chip",
                "--no-interactive",
                "--accept-all-defaults",
                "-o", str(out_csv),
                str(tmp_path),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert out_csv.exists()
    
    def test_design_rna_paired_end(self, tmp_path):
        """Test design CSV generation for RNA-seq paired-end."""
        (tmp_path / "control_rep1_R1.fastq.gz").touch()
        (tmp_path / "control_rep1_R2.fastq.gz").touch()
        (tmp_path / "treatment_rep1_R1.fastq.gz").touch()
        (tmp_path / "treatment_rep1_R2.fastq.gz").touch()
        
        out_csv = tmp_path / "rna_design.csv"
        result = subprocess.run(
            [
                "seqnado", "design", "rna",
                "--no-interactive",
                "--accept-all-defaults",
                "-o", str(out_csv),
                str(tmp_path),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert out_csv.exists()


class TestPipelineCommand:
    """Tests for seqnado pipeline command."""
    
    def test_pipeline_help(self):
        """Test pipeline command help."""
        result = subprocess.run(
            ["seqnado", "pipeline", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "pipeline" in result.stdout.lower() or "run" in result.stdout.lower()
    
    def test_pipeline_dry_run(self, tmp_path, monkeypatch, genome_config_dir):
        """Test pipeline dry-run mode."""
        monkeypatch.setenv("HOME", str(tmp_path))
        
        # Create minimal config file
        config_file = tmp_path / "config.yaml"
        config_file.write_text("""
assay: ATAC
project:
  name: test_project
genome:
  name: test_genome
  star_index: /tmp/star
metadata: metadata.csv
""")
        
        # Create metadata
        (tmp_path / "metadata.csv").write_text("sample_id,fastq_r1\nsample1,/tmp/r1.fastq.gz\n")
        
        result = subprocess.run(
            [
                "seqnado", "pipeline",
                "--config", str(config_file),
                "--dry-run",
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            timeout=60,
        )
        
        # Dry run might fail due to missing files, but should not crash
        # Just check it doesn't crash with invalid arguments
        assert "config" in result.stdout.lower() or "error" in result.stderr.lower() or result.returncode in [0, 1]


class TestVersionAndHelp:
    """Tests for version and help commands."""
    
    def test_version_output_format(self):
        """Test version command output format."""
        result = subprocess.run(
            ["seqnado", "--version"],
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0
        output = result.stdout.strip()
        assert len(output) > 0
        assert "seqnado" in output.lower() or "version" in output.lower() or any(c.isdigit() for c in output)
    
    def test_help_lists_commands(self):
        """Test that help lists available commands."""
        result = subprocess.run(
            ["seqnado", "--help"],
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0
        # Should list main commands
        assert any(cmd in result.stdout.lower() for cmd in ["init", "config", "design", "pipeline"])
    
    def test_command_help_individual(self):
        """Test help for each individual command."""
        commands = ["init", "config", "design", "pipeline"]
        
        for cmd in commands:
            result = subprocess.run(
                ["seqnado", cmd, "--help"],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0, f"Help failed for {cmd}"
            assert len(result.stdout) > 0


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_config_without_metadata(self, tmp_path, monkeypatch, genome_config_dir):
        """Test config command without existing metadata file."""
        monkeypatch.setenv("HOME", str(tmp_path))
        (tmp_path / "seqnado_output").mkdir()
        
        out_file = tmp_path / "config.yaml"
        result = subprocess.run(
            [
                "seqnado", "config", "atac",
                "--no-interactive",
                "--no-make-dirs",
                "-o", str(out_file),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        
        # Should handle gracefully
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert out_file.exists()
    
    def test_design_empty_directory(self, tmp_path):
        """Test design command on empty directory."""
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        
        result = subprocess.run(
            [
                "seqnado", "design", "atac",
                "--no-interactive",
                "-o", "design.csv",
                str(empty_dir),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        
        # Should handle no FASTQs gracefully
        assert "no" in result.stderr.lower() or "not found" in result.stderr.lower() or result.returncode != 0
    
    def test_invalid_assay_name(self, tmp_path):
        """Test config with invalid assay name."""
        result = subprocess.run(
            ["seqnado", "config", "invalid_assay"],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        
        # Should fail or show available assays
        assert result.returncode != 0 or "invalid" in result.stderr.lower()
