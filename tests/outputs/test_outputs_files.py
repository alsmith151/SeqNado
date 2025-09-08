from pathlib import Path

import pytest

from seqnado import Assay, PileupMethod, DataScalingTechnique, PeakCallingMethod, MethylationMethod
from seqnado.outputs.files import (
    BigWigFiles,
    PeakCallingFiles,
    HeatmapFiles,
    HubFiles,
    PlotFiles,
    SpikeInFiles,
    SNPFilesRaw,
    SNPFilesAnnotated,
    MethylationFiles,
    BigBedFiles,
)


def test_bigwigfiles_rna_stranded_paths():
    bw = BigWigFiles(
        assay=Assay.RNA,
        names=["s1"],
        pileup_methods=[PileupMethod.DEEPTOOLS],
        scale_methods=[DataScalingTechnique.UNSCALED, DataScalingTechnique.CPM],
    )
    files = bw.files
    # 2 scales x 2 strands = 4 files
    assert len(files) == 4
    assert any(f.endswith("/deeptools/unscaled/s1_plus.bigWig") for f in files)
    assert any(f.endswith("/deeptools/unscaled/s1_minus.bigWig") for f in files)
    assert any(f.endswith("/deeptools/cpm/s1_plus.bigWig") for f in files)
    assert any(f.endswith("/deeptools/cpm/s1_minus.bigWig") for f in files)


def test_bigwigfiles_skip_incompatible():
    bw = BigWigFiles(
        assay=Assay.ATAC,
        names=["s1"],
        pileup_methods=[PileupMethod.DEEPTOOLS, PileupMethod.HOMER],
        scale_methods=[DataScalingTechnique.UNSCALED, DataScalingTechnique.SPIKEIN],
    )
    files = bw.files
    # HOMER + SPIKEIN is incompatible, so expect 3 files total
    assert len(files) == 3
    assert any(f.endswith("/deeptools/unscaled/s1.bigWig") for f in files)
    assert any(f.endswith("/deeptools/spikein/s1.bigWig") for f in files)
    assert any(f.endswith("/homer/unscaled/s1.bigWig") for f in files)
    assert not any("/homer/spikein/" in f for f in files)


def test_peakcallingfiles_paths():
    peaks = PeakCallingFiles(
        assay=Assay.ATAC,
        names=["s1", "s2"],
        peak_calling_method=[PeakCallingMethod.MACS2, PeakCallingMethod.HOMER],
    )
    files = peaks.files
    assert "seqnado_output/peaks/macs2/s1.bed" in files
    assert "seqnado_output/peaks/macs2/s2.bed" in files
    assert "seqnado_output/peaks/homer/s1.bed" in files
    assert "seqnado_output/peaks/homer/s2.bed" in files


def test_heatmapfiles_invalid_assay_raises():
    with pytest.raises(ValueError):
        HeatmapFiles(assay=Assay.SNP)


def test_hubfiles_paths(tmp_path: Path):
    hub_dir = tmp_path / "hub"
    name = "myhub"
    hf = HubFiles(hub_dir=hub_dir, hub_name=name)
    assert hf.files == [str(hub_dir / f"{name}.hub.txt")]


def test_plotfiles_reads_coordinates(tmp_path: Path):
    # Create a minimal BED-like coordinates file
    coords = tmp_path / "coords.bed"
    coords.write_text("chr1\t100\t200\tregionA\nchr2\t300\t400\tregionB\n")

    pf = PlotFiles(coordinates=coords, file_format="pdf")
    files = pf.files
    assert len(files) == 2
    # names taken from 4th column when present
    assert "seqnado_output/genome_browser_plots/regionA.pdf" in [str(f) for f in files]
    assert "seqnado_output/genome_browser_plots/regionB.pdf" in [str(f) for f in files]


def test_spikeinfiles_invalid_assay_raises():
    with pytest.raises(ValueError):
        SpikeInFiles(assay=Assay.SNP, names=["s1"])  # SNP not in AssaysWithSpikein


def test_snp_files_paths():
    raw = SNPFilesRaw(assay=Assay.SNP, names=["a", "b"]).files
    anno = SNPFilesAnnotated(assay=Assay.SNP, names=["a", "b"]).files
    assert raw == [
        "seqnado_output/variant/a.vcf.gz",
        "seqnado_output/variant/b.vcf.gz",
    ]
    assert anno == [
        "seqnado_output/variant/a.anno.vcf.gz",
        "seqnado_output/variant/b.anno.vcf.gz",
    ]


def test_methylationfiles_paths():
    mf = MethylationFiles(
        assay=Assay.METH,
        names=["s1"],
        genomes=["hg38", "dm6"],
        method=MethylationMethod.TAPS,
    )
    files = mf.files
    # split BAMs
    assert "seqnado_output/aligned/spikein/s1_hg38.bam" in files
    assert "seqnado_output/aligned/spikein/s1_dm6.bam" in files
    # methyldackel outputs include method in name
    assert "seqnado_output/methylation/methyldackel/s1_hg38_tapsCpG.bedGraph" in files
    assert "seqnado_output/methylation/methyldackel/s1_dm6_tapsCpG.bedGraph" in files
    # methylation conversion + bias files
    assert "seqnado_output/methylation/methylation_conversion.tsv" in files
    assert "seqnado_output/methylation/methyldackel/bias/s1_hg38.txt" in files
    assert "seqnado_output/methylation/methyldackel/bias/s1_dm6.txt" in files


def test_bigbedfiles_from_bed_list(tmp_path: Path):
    bed1 = tmp_path / "a.bed"
    bed1.write_text("chr1\t1\t10\n")
    txt = tmp_path / "b.txt"
    txt.write_text("x\n")

    bb = BigBedFiles(bed_files=[bed1, txt])
    files = bb.files
    assert files == [str(bed1.with_suffix(".bb"))]
