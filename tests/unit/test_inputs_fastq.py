from __future__ import annotations

import csv
from pathlib import Path

from seqnado import Assay
from seqnado.inputs.core import clean_sample_name
from seqnado.inputs.fastq import FastqCollection, FastqFile, FastqFileIP, FastqSet
from seqnado.inputs.helpers import get_sample_collection


def _write_fastq(tmp: Path, name: str) -> Path:
    p = tmp / name
    p.write_text("@r\nN\n+\n#\n")
    return p


def test_fastqfile_parses_read_and_sample(tmp_path: Path):
    f = _write_fastq(tmp_path, "chip_A_S1_L001_R1_001.fastq.gz")
    fq = FastqFile(path=f)
    assert fq.read_number == 1
    assert fq.is_paired is True
    assert clean_sample_name("chip_A_S1_L001_R1_001") == fq.sample_base


def test_fastqfileip_predicts_ip_and_control(tmp_path: Path):
    f = _write_fastq(tmp_path, "sampleX_IGG_L001_R1_001.fastq.gz")
    fq = FastqFileIP(path=f)
    assert fq.ip == "IGG"
    assert fq.is_control is True
    assert fq.sample_base_without_ip == "sampleX"


def test_fastqset_single_and_paired(tmp_path: Path):
    r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
    fs_single = FastqSet(r1=r1)
    assert fs_single.is_paired is False
    assert fs_single.file_paths == [r1.path]

    r2 = FastqFile(path=_write_fastq(tmp_path, "s1_R2.fastq.gz"))
    fs_pair = FastqSet(r1=r1, r2=r2)
    assert fs_pair.is_paired is True
    assert set(fs_pair.file_paths) == {r1.path, r2.path}


def test_fastqcollection_symlink(tmp_path: Path):
    # Create simple paired sample
    r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
    r2 = FastqFile(path=_write_fastq(tmp_path, "s1_R2.fastq.gz"))
    fc = FastqCollection(
        assay=Assay.RNA,
        metadata=[],
        fastq_sets=[FastqSet(sample_id="s1", r1=r1, r2=r2)],
    )

    out = tmp_path / "links"
    fc.symlink_fastq_files(out)
    assert (out / "s1_1.fastq.gz").is_symlink()
    assert (out / "s1_2.fastq.gz").is_symlink()


def test_get_sample_collection_from_csv(tmp_path: Path):
    # Realistic metadata CSV with R1/R2 columns triggers FastqCollection
    csv_path = tmp_path / "metadata.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "r1", "r2"])
        w.writerow(
            [
                "s1",
                str(_write_fastq(tmp_path, "s1_R1.fastq.gz")),
                str(_write_fastq(tmp_path, "s1_R2.fastq.gz")),
            ]
        )

    sc = get_sample_collection(Assay.RNA, csv_path)
    # The returned type implements CollectionLike; check basic behavior
    paths = sc.get_file_paths("fastq")
    assert len(paths) == 2
    assert all(Path(p).exists() for p in paths)
