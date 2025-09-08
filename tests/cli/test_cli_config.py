from __future__ import annotations

from pathlib import Path
from click.testing import CliRunner
import json
import builtins


def _write_minimal_genome_config(tmp_root: Path) -> Path:
    """Create a minimal genome_config.json that passes validation for RNA (STAR index).
    Includes a small fasta file to satisfy required key and existence checks.
    """
    cfg_dir = tmp_root / ".config" / "seqnado"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    star_dir = tmp_root / "star_index"
    star_dir.mkdir(parents=True, exist_ok=True)
    fasta = tmp_root / "genome.fa"
    fasta.write_text(">chr1\n" "NNNNNNNNNNNNNNNNNNNN\n")

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
    # Arrange workspace: minimal genome config, metadata file, and template
    _write_minimal_genome_config(tmp_path)
    (tmp_path / "metadata.csv").write_text("sample,fastq\n")
    # No need to copy template; cli now uses packaged resource

    # Ensure CLI reads genome config from our tmp root
    monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))
    # Fake package version lookup
    monkeypatch.setenv("PYTHONHASHSEED", "0")  # reduce flakiness
    def fake_version(_):
        return "0.0.0-test"
    import importlib.metadata as _md
    monkeypatch.setattr(_md, "version", fake_version)

    # Avoid file lookup inside get_tool_options by returning a minimal string
    import seqnado.config.user_input as ui
    monkeypatch.setattr(ui, "get_tool_options", lambda assay: "options: {}\n")

    # Simulate pressing Enter for every prompt to accept defaults
    monkeypatch.setattr(builtins, "input", lambda _: "")

    # Change to temp working directory
    monkeypatch.chdir(tmp_path)

    # Act: invoke CLI to generate config for RNA, without creating extra directories
    from seqnado.cli import cli_config
    runner = CliRunner()
    res = runner.invoke(cli_config, ["rna", "--dont-make-directories"])  # type: ignore[arg-type]

    # Assert
    assert res.exit_code == 0, res.output
    out_file = tmp_path / "config_RNA.yaml"
    assert out_file.exists(), "Config file should be created"
    text = out_file.read_text()
    # Light sanity checks on content
    assert "Project" in text or "project" in text
    assert "metadata.csv" in text
    assert "assay_config" in text


def test_cli_config_rna_all_options_flag(monkeypatch, tmp_path: Path):
    # Arrange: environment and inputs as above
    _write_minimal_genome_config(tmp_path)
    (tmp_path / "metadata.csv").write_text("sample,fastq\n")
    # No need to copy template; cli now uses packaged resource

    monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))
    def fake_version(_):
        return "0.0.0-test"
    import importlib.metadata as _md
    monkeypatch.setattr(_md, "version", fake_version)

    import seqnado.config.user_input as ui
    monkeypatch.setattr(builtins, "input", lambda _: "")
    monkeypatch.chdir(tmp_path)

    from seqnado.cli import cli_config
    runner = CliRunner()
    res = runner.invoke(cli_config, ["rna", "--dont-make-directories", "--all-options"])  # type: ignore[arg-type]

    assert res.exit_code == 0, res.output
    out_file = tmp_path / "config_RNA.yaml"
    assert out_file.exists()
