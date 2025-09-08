from pathlib import Path

from seqnado.inputs import Assay, BamCollection


def _touch(p: Path) -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("")
    return p


def test_bamcollection_from_files_and_dataframe(tmp_path: Path):
    f1 = _touch(tmp_path / "a.bam")
    f2 = _touch(tmp_path / "b.bam")

    def md_builder(sample: str):
        from seqnado.inputs import Metadata

        return Metadata(consensus_group="g1")

    coll = BamCollection.from_files(
        assay=Assay.ATAC, files=[f1, f2], metadata=md_builder
    )
    assert coll.primary_file_type == "bam"
    assert sorted(coll.sample_ids) == ["a", "b"]
    assert set(coll.get_file_paths()) == {f1, f2}
    df = coll.to_dataframe()
    # Expect rows for both samples and metadata columns present
    assert set(df["sample_id"]) == {"a", "b"}
    assert "consensus_group" in df.columns and "assay" in df.columns


def test_bamcollection_get_file_paths_kinds_and_symlinks(tmp_path: Path):
    f1 = _touch(tmp_path / "s1.bam")
    coll = BamCollection.from_files(assay=Assay.ATAC, files=[f1])
    # kind filter
    assert coll.get_file_paths("bam") == [f1]
    try:
        coll.get_file_paths("fastq")
        raised = False
    except ValueError:
        raised = True
    assert raised
    # symlink
    out = tmp_path / "links"
    coll.symlink_bam_files(out)
    assert (out / "s1.bam").exists()
