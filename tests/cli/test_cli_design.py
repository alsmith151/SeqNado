from pathlib import Path
import subprocess


def _write_fastq(tmp: Path, name: str) -> Path:
    p = tmp / name
    p.touch()
    return p


def test_cli_design_generates_csv(tmp_path: Path, monkeypatch):
    # create a couple fastqs
    r1 = _write_fastq(tmp_path, "chipA_H3K27ac_L001_R1_001.fastq.gz")
    r2 = _write_fastq(tmp_path, "chipA_H3K27ac_L001_R2_001.fastq.gz")

    # Invoke CLI via subprocess with explicit file paths
    result = subprocess.run(
        [
            "seqnado",
            "design",
            "chip",
            "-o",
            "metadata.csv",
            "--no-interactive",
            "--accept-all-defaults",
            str(r1),
            str(r2),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
    out = tmp_path / "metadata.csv"
    assert out.exists(), "metadata.csv should be created"
    text = out.read_text()
    # basic sanity checks on columns
    assert "sample" in text.lower()  # Could be sample_id or sample_name
    assert "r1" in text.lower()