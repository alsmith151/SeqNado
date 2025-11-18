import subprocess


def test_seqnado_help_exits_zero():
    res = subprocess.run(["seqnado", "--help"], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert "SeqNado" in res.stdout or "Usage" in res.stdout


def test_seqnado_pipeline_help_exits_zero():
    res = subprocess.run(["seqnado", "pipeline", "--help"], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert "pipeline" in res.stdout.lower() or "Usage" in res.stdout


def test_seqnado_version_exits_zero():
    res = subprocess.run(["seqnado", "--version"], capture_output=True, text=True)
    assert res.returncode == 0
    assert res.stdout.strip() != ""