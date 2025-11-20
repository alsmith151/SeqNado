from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess

import pytest

from seqnado.config.core import SeqnadoConfig


def _tests_data_genome(repo_root: Path) -> Path:
    return repo_root / "tests" / "data" / "genome"


def _find_star_or_bt2(genome_root: Path) -> tuple[str, Path | None, Path | None]:
    """Return (assay, star_dir, bt2_prefix) using tests/data/genome if available.

    Prioritize STAR (assay='rna') if a STAR* directory exists; else look for any .bt2 file
    and compute the Bowtie2 prefix path by stripping typical suffixes.
    If nothing found, return (None, None, None).
    """
    if not genome_root.exists():
        return (None, None, None)

    # STAR index directory candidate
    for p in genome_root.iterdir():
        if p.is_dir() and p.name.lower().startswith("star"):
            return ("rna", p, None)

    # Bowtie2 index: look for any *.bt2 under subdirectories and infer prefix
    import re
    for bt2 in genome_root.rglob("*.bt2*"):
        name = bt2.name
        # remove .1.bt2, .2.bt2, .rev.1.bt2(l), etc.
        prefix_name = re.sub(r"\.(?:[1-4]|rev\.[12])\.bt2l?$", "", name)
        prefix_path = bt2.parent / prefix_name
        return ("chip", None, prefix_path)

    return (None, None, None)


def _write_genome_config(home: Path, assay: str, star_dir: Path | None, bt2_prefix: Path | None) -> Path:
    cfg_dir = home / ".config" / "seqnado"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    # Use a real genome key so organism prediction works in model
    key = "hg38"

    payload: dict[str, dict] = {key: {}}
    if assay == "rna" and star_dir is not None:
        payload[key]["star_index"] = str(star_dir)
    elif assay == "chip" and bt2_prefix is not None:
        # Bowtie2 expects a prefix without suffix; validation will check matching files exist
        payload[key]["bt2_index"] = str(bt2_prefix)
    cfg = cfg_dir / "genome_config.json"
    cfg.write_text(json.dumps(payload, indent=2))
    return cfg


@pytest.mark.integration
def test_cli_config_roundtrip_with_tests_data(tmp_path: Path):
    repo_root = Path(__file__).resolve().parents[2]
    genome_root = _tests_data_genome(repo_root)
    assay, star_dir, bt2_dir = _find_star_or_bt2(genome_root)

    if assay is None:
        pytest.skip("No STAR or Bowtie2 test data found under tests/data/genome")

    # Write genome config under isolated HOME/SEQNADO_CONFIG
    os.environ["HOME"] = str(tmp_path)
    os.environ["SEQNADO_CONFIG"] = str(tmp_path)
    _write_genome_config(tmp_path, assay, star_dir, bt2_dir)

    # Create seqnado_output dir for UCSCHubConfig validation
    (tmp_path / "seqnado_output").mkdir()

    # Generate config via CLI (non-interactive) into a known file
    out = tmp_path / f"config_{assay}.yaml"
    res = subprocess.run(
        [
            "seqnado",
            "config",
            assay,
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

    if res.returncode != 0:
        print("STDOUT:\n", res.stdout)
        print("STDERR:\n", res.stderr)
    assert res.returncode == 0
    assert out.exists()

    # Load config into model
    cfg = SeqnadoConfig.from_yaml(out)
    assert cfg.assay.value.lower() == assay
    # Roundtrip: dump to YAML-equivalent dict and reload
    import yaml
    dump = tmp_path / "roundtrip.yaml"
    with dump.open("w") as f:
        yaml.safe_dump(cfg.model_dump(mode="json"), f)
    cfg2 = SeqnadoConfig.from_yaml(dump)
    assert cfg2.assay == cfg.assay
    assert cfg2.project.name == cfg.project.name