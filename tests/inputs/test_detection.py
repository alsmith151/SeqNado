import pytest

from seqnado.inputs.interfaces import detect_file_type


@pytest.mark.parametrize(
    "paths,expected",
    [
        (["a.fastq", "b.fastq.gz", "c.FQ"], "fastq"),
        (["x.bam"], "bam"),
        (["x.bigWig", "y.bw"], "bigwig"),
        (["unknown.txt", "notes.md"], None),
        # Majority rule
        (["a.fastq", "b.bam", "c.fastq.gz"], "fastq"),
        (["a.bam", "b.bam", "c.bw"], "bam"),
        # Tie -> None
        (["a.fastq", "b.bam"], None),
        (["a.bw", "b.bigWig"], "bigwig"),
    ],
)
def test_detect_file_type_majority(paths, expected):
    assert detect_file_type(paths) == expected
