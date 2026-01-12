from __future__ import annotations

import json
import subprocess
from pathlib import Path


def _write_minimal_genome_config(tmp_root: Path) -> Path:
    """Create a minimal genome_config.json that passes validation for RNA (STAR index).
    Includes a small fasta file to satisfy required key and existence checks.
    """
    cfg_dir = tmp_root / ".config" / "seqnado"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    star_dir = tmp_root / "star_index"
    star_dir.mkdir(parents=True, exist_ok=True)
    fasta = tmp_root / "genome.fa"
    fasta.write_text(">chr1\nNNNNNNNNNNNNNNNNNNNN\n")

    genome_cfg = {
        "test": {
            # Only include the required STAR index for RNA; omit optional paths to avoid path validation
            "star_index": str(star_dir),
            "fasta": str(fasta),
        }
    }

    cfg_path = cfg_dir / "genome_config.json"
    cfg_path.write_text(json.dumps(genome_cfg, indent=2))
    return cfg_path


def test_cli_config_rna_creates_config_file(monkeypatch, tmp_path: Path):
    # Arrange workspace: minimal genome config, metadata file
    _write_minimal_genome_config(tmp_path)
    (tmp_path / "metadata.csv").write_text("sample_id,fastq\n")

    # Ensure CLI reads genome config from our tmp root
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))

    # Create seqnado_output dir for UCSCHubConfig validation
    (tmp_path / "seqnado_output").mkdir()

    # Act: invoke CLI via subprocess to generate config for RNA
    out_file = tmp_path / "config_rna.yaml"
    result = subprocess.run(
        [
            "seqnado",
            "config",
            "rna",
            "--no-interactive",
            "--no-make-dirs",
            "--render-options",
            "-o",
            str(out_file),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    # Assert
    assert result.returncode == 0, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
    assert out_file.exists(), "Config file should be created"
    text = out_file.read_text()
    # Light sanity checks on content
    assert "project" in text.lower()
    assert (
        "metadata" in text
    )  # Check for metadata field (could be metadata.csv or metadata_rna.csv)
    assert "assay_config" in text


def test_cli_config_rna_all_options_flag(monkeypatch, tmp_path: Path):
    # Arrange: environment and inputs
    _write_minimal_genome_config(tmp_path)
    (tmp_path / "metadata.csv").write_text("sample_id,fastq\n")

    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))

    # Create seqnado_output dir for UCSCHubConfig validation
    (tmp_path / "seqnado_output").mkdir()

    # Act: invoke CLI with --render-options flag
    out_file = tmp_path / "config_rna.yaml"
    result = subprocess.run(
        [
            "seqnado",
            "config",
            "rna",
            "--no-interactive",
            "--no-make-dirs",
            "--render-options",
            "-o",
            str(out_file),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    # Assert
    assert result.returncode == 0, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
    assert out_file.exists()
    # When --render-options is used, the config should include detailed comments/options
    text = out_file.read_text()
    assert len(text) > 100  # Should have substantial content with options
