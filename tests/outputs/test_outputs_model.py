from seqnado.outputs.core import SeqnadoOutputFiles
from seqnado import DataScalingTechnique, PileupMethod


def test_seqnado_output_files_filters():
    files = [
        "seqnado_output/bw/deeptools/unscaled/a.bigWig",
        "seqnado_output/bw/deeptools/cpm/a.bigWig",
        "seqnado_output/bw/homer/unscaled/b.bigWig",
        "seqnado_output/peaks/macs2/a.bed",
        "seqnado_output/peaks/macs2/b.bed",
        "seqnado_output/peaks/merged/consensus.bed",
        "seqnado_output/bigbeds/x.bigBed",
        "seqnado_output/heatmap/heatmap.pdf",
        "seqnado_output/genome_browser_plots/region.pdf",
    ]

    out = SeqnadoOutputFiles(files=files, sample_names=["a", "b"]) 

    # bigwig files
    assert len(out.bigwig_files) == 3

    # select subtype
    sel = out.select_bigwig_subtype(method=PileupMethod.DEEPTOOLS, scale=DataScalingTechnique.CPM)
    assert sel == ["seqnado_output/bw/deeptools/cpm/a.bigWig"]

    # peaks
    assert set(out.peak_files) == {"seqnado_output/peaks/macs2/a.bed", "seqnado_output/peaks/macs2/b.bed", "seqnado_output/peaks/merged/consensus.bed"}
    assert out.has_consensus_peaks is True

    # bigbed
    assert out.bigbed_files == ["seqnado_output/bigbeds/x.bigBed"]

    # heatmaps
    assert out.heatmap_files == ["seqnado_output/heatmap/heatmap.pdf"]

    # genome browser plots
    assert out.genome_browser_plots == ["seqnado_output/genome_browser_plots/region.pdf"]
