from pathlib import Path
import pytest

from seqnado.config import STARIndex


def test_star_index_requires_existing_dir(tmp_path: Path):
    with pytest.raises(ValueError):
        STARIndex(prefix=tmp_path / "missing")
    # Valid when dir exists
    d = tmp_path / "ok"
    d.mkdir()
    s = STARIndex(prefix=d)
    assert isinstance(s, STARIndex)
    assert s.type == "STAR"


def test_bowtie_index_validates_and_lists_files(tmp_path: Path):
    # Create bowtie2 index file set for prefix
    prefix = tmp_path / "bt2" / "genome"
    prefix.parent.mkdir(parents=True)
    suffixes = [
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    ]
    for sfx in suffixes:
        (Path(str(prefix) + sfx)).write_text("")

    from seqnado.config import BowtieIndex

    idx = BowtieIndex(prefix=str(prefix))
    files = idx.files
    assert len(files) == 6 and all(f.exists() for f in files)
