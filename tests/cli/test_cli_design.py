from pathlib import Path
from click.testing import CliRunner


def _write_fastq(tmp: Path, name: str) -> Path:
    p = tmp / name
    p.write_text("@r\nN\n+\n#\n")
    return p


def test_cli_design_generates_csv(tmp_path: Path, monkeypatch):
    # create a couple fastqs in cwd so cli can auto-discover if files arg omitted
    _write_fastq(tmp_path, "chipA_H3K27ac_L001_R1_001.fastq.gz")
    _write_fastq(tmp_path, "chipA_H3K27ac_L001_R2_001.fastq.gz")
    monkeypatch.chdir(tmp_path)

    # import the command and call via CliRunner
    from seqnado.cli import cli_design
    r = CliRunner().invoke(cli_design, ["chip", "--merge", "-o", "design.csv"])  # type: ignore[arg-type]

    assert r.exit_code == 0, r.output
    out = tmp_path / "design.csv"
    text = out.read_text()
    # basic sanity checks on columns
    assert "sample_name" in text.lower()
    assert "r1" in text.lower()
    assert "merge" in text.lower()