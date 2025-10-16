from pathlib import Path

import pytest

from seqnado.config.configs import GenomeConfig, BowtieIndex, STARIndex, UserFriendlyError


def test_bowtie_index_files_validation(tmp_path: Path):
    # Create fake bowtie2 index files with typical suffixes
    idx_prefix = tmp_path / "idx/genome"
    idx_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in [".1.bt2", ".2.bt2", ".rev.1.bt2"]:
        (Path(str(idx_prefix) + suffix)).write_text("")

    cfg = GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(idx_prefix)))
    assert isinstance(cfg.index, BowtieIndex)
    assert cfg.index.files, "Expected to find bowtie index files"


def test_star_index_validation(tmp_path: Path):
    idx_dir = tmp_path / "star"
    idx_dir.mkdir(parents=True)
    (idx_dir / "SA").write_text("")  # dummy content
    cfg = GenomeConfig(name="hg38", index=STARIndex(prefix=idx_dir))
    assert isinstance(cfg.index, STARIndex)


def test_missing_optional_paths_raise_userfriendlyerror(tmp_path: Path):
    # fasta points to a missing file -> should raise our friendly error
    with pytest.raises(UserFriendlyError):
        GenomeConfig(name="hg38", index=STARIndex(prefix=tmp_path), fasta=tmp_path / "missing.fa")
