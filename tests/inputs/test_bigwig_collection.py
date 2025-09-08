from pathlib import Path

from seqnado.inputs import Assay, BigWigCollection, ensure_collection


def _touch(p: Path) -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("")
    return p


def test_bigwigcollection_rna_strand_logic_and_dataframe(tmp_path: Path):
    s1_plus = _touch(tmp_path / "s1_plus.bw")
    s1_minus = _touch(tmp_path / "s1_minus.bw")
    s2_plus = _touch(tmp_path / "s2_plus.bigWig")
    s2_minus = _touch(tmp_path / "s2_minus.bigWig")

    coll = BigWigCollection.from_files(
        assay=Assay.RNA,
        files=[s1_minus, s2_plus, s1_plus, s2_minus],
    )
    # sample ids derive from stem without strand suffix
    assert sorted(coll.sample_ids) == ["s1", "s1", "s2", "s2"]
    # tracks map to lists ordered plus before minus
    tracks = coll.tracks
    assert [p.stem for p in tracks["s1"]] == ["s1_plus", "s1_minus"]
    assert [p.stem for p in tracks["s2"]] == ["s2_plus", "s2_minus"]
    # dataframe has plus/minus columns
    df = coll.to_dataframe()
    assert set(df["sample_id"]) == {"s1", "s2"}
    assert {"bigwig_plus", "bigwig_minus"}.issubset(df.columns)


def test_bigwigcollection_non_rna_single_track_and_symlinks(tmp_path: Path):
    s1 = _touch(tmp_path / "s1.bw")
    s2 = _touch(tmp_path / "s2.bigWig")
    coll = BigWigCollection.from_files(assay=Assay.ATAC, files=[s1, s2])
    assert coll.primary_file_type == "bigwig"
    assert sorted(coll.sample_names) == ["s1", "s2"]
    df = coll.to_dataframe()
    assert "bigwig" in df.columns and set(df["sample_id"]) == {"s1", "s2"}
    out = tmp_path / "bw_links"
    coll.symlink_bigwig_files(out)
    assert (out / s1.name).exists() and (out / s2.name).exists()


def test_ensure_collection_accepts_bigwig(tmp_path: Path):
    bw = BigWigCollection.from_files(assay=Assay.RNA, files=[_touch(tmp_path / "s_plus.bw")])
    assert ensure_collection(bw) is bw
