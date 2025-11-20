import os
from pathlib import Path
import json
import subprocess


def _write_minimal_star_genome_config(root: Path) -> Path:
    cfg_dir = root / ".config" / "seqnado"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    star_dir = root / "star_index"
    star_dir.mkdir(parents=True, exist_ok=True)

    genome_cfg = {
        "test": {
            "star_index": str(star_dir),
            # keep other fields unset to avoid validation of non-existing paths
        }
    }
    cfg_path = cfg_dir / "genome_config.json"
    cfg_path.write_text(json.dumps(genome_cfg, indent=2))
    return cfg_path


def test_cli_config_rna_non_interactive_no_makedirs(tmp_path: Path, monkeypatch):
    # Arrange: minimal STAR genome config in temp HOME/SEQNADO_CONFIG
    _write_minimal_star_genome_config(tmp_path)
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))

    # Create seqnado_output dir so UCSCHubConfig validator doesn't fail
    (tmp_path / "seqnado_output").mkdir()

    out = tmp_path / "config_RNA.yaml"
    res = subprocess.run(
        [
            "seqnado",
            "config",
            "rna",
            "--no-interactive",
            "--no-make-dirs",
            "--render-options",
            "-o",
            str(out),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, f"stderr=\n{res.stderr}\nstdout=\n{res.stdout}"
    assert out.exists(), "Rendered config file should exist"
    text = out.read_text()
    assert "assay" in text.lower()
    assert "genome" in text.lower() or "project" in text.lower()


def test_cli_config_rna_non_interactive_creates_project_dir(tmp_path: Path, monkeypatch):
    _write_minimal_star_genome_config(tmp_path)
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))

    # Create seqnado_output dir for UCSCHubConfig validation
    (tmp_path / "seqnado_output").mkdir()

    res = subprocess.run(
        [
            "seqnado",
            "config",
            "rna",
            "--no-interactive",
            "--render-options",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, f"stderr=\n{res.stderr}\nstdout=\n{res.stdout}"

    # Expect a directory named YYYY-MM-DD_<project>/config_rna.yaml
    candidates = list(tmp_path.glob("*/config_rna.yaml"))
    assert candidates, "Expected config_rna.yaml inside a dated project directory"