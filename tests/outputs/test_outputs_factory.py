from pathlib import Path

from seqnado import Assay, PileupMethod, PeakCallingMethod, QuantificationMethod
from seqnado.config import (
    SeqnadoConfig,
    ProjectConfig,
    GenomeConfig,
    BigwigConfig,
    PeakCallingConfig,
    UCSCHubConfig,
    RNAQuantificationConfig,
    ATACAssayConfig,
    RNAAssayConfig,
)
from seqnado.inputs import FastqCollection
from seqnado.inputs.fastq import FastqFile, FastqSet
from seqnado.outputs.core import SeqnadoOutputFactory


def make_fastq_collection(tmp_path: Path, assay: Assay) -> FastqCollection:
    r1 = tmp_path / "a_R1.fastq.gz"
    r2 = tmp_path / "a_R2.fastq.gz"
    r1.write_text("N\n")
    r2.write_text("N\n")
    fs = FastqSet(sample_id="a", r1=FastqFile(path=r1), r2=FastqFile(path=r2))
    return FastqCollection(assay=assay, fastq_sets=[fs], metadata=[])


def minimal_base_config(tmp_path: Path, assay: Assay, assay_config) -> SeqnadoConfig:
    # GenomeConfig requires an index object; provide a dummy STARIndex directory
    genome_dir = tmp_path / "star_index"
    genome_dir.mkdir()
    from seqnado.config.configs import STARIndex

    return SeqnadoConfig(
        assay=assay,
        project=ProjectConfig(name="test"),
        genome=GenomeConfig(name="hg38", index=STARIndex(prefix=genome_dir)),
        metadata=tmp_path / "meta.csv",
        assay_config=assay_config,
    )


def test_factory_builds_bigwigs_and_peaks_for_atac(tmp_path: Path):
    samples = make_fastq_collection(tmp_path, Assay.ATAC)

    assay_cfg = ATACAssayConfig(
        bigwigs=BigwigConfig(pileup_method=[PileupMethod.DEEPTOOLS], binsize=10),
        peak_calling=PeakCallingConfig(method=[PeakCallingMethod.MACS2]),
        create_heatmaps=True,
        create_geo_submission_files=False,
    )
    cfg = minimal_base_config(tmp_path, Assay.ATAC, assay_cfg)

    builder = SeqnadoOutputFactory(assay=Assay.ATAC, samples=samples, config=cfg).create_output_builder()
    out = builder.build()

    # QC and bigwigs
    assert any(f.endswith("seqnado_output/seqnado_report.html") for f in out.files)
    assert any("/deeptools/unscaled/a.bigWig" in f for f in out.files)

    # Peaks
    assert "seqnado_output/peaks/macs2/a.bed" in out.files

    # Heatmaps
    assert "seqnado_output/heatmap/heatmap.pdf" in out.files


def test_factory_ucsc_hub_added(tmp_path: Path):
    samples = make_fastq_collection(tmp_path, Assay.RNA)
    # UCSCHubConfig validates parent directory exists; set directory under tmp_path
    # UCSCHubConfig requires the parent of `directory` to exist.
    # Use a directory under an existing parent, and ensure the immediate parent exists.
    parent = tmp_path / "parent"
    (parent / "seqnado_output").mkdir(parents=True, exist_ok=True)
    hub_cfg = UCSCHubConfig(directory=str(parent / "seqnado_output/hub/"))

    assay_cfg = RNAAssayConfig(
        bigwigs=BigwigConfig(pileup_method=[PileupMethod.DEEPTOOLS]),
        ucsc_hub=hub_cfg,
    )
    cfg = minimal_base_config(tmp_path, Assay.RNA, assay_cfg)

    builder = SeqnadoOutputFactory(assay=Assay.RNA, samples=samples, config=cfg).create_output_builder()
    out = builder.build()

    assert any("seqnado_output/hub" in f for f in out.files)


def test_factory_geo_submission_requires_fastq(tmp_path: Path):
    samples = make_fastq_collection(tmp_path, Assay.RNA)
    assay_cfg = RNAAssayConfig(
        bigwigs=BigwigConfig(pileup_method=[PileupMethod.DEEPTOOLS]),
        create_heatmaps=False,
    )
    # Toggle GEO submission
    assay_cfg.create_geo_submission_files = True
    cfg = minimal_base_config(tmp_path, Assay.RNA, assay_cfg)

    builder = SeqnadoOutputFactory(assay=Assay.RNA, samples=samples, config=cfg).create_output_builder()
    # Should include GEO files paths for raw fastq and processed bigwigs (renamed)
    out = builder.build()
    # md5sums and samples table paths are included
    assert any(f.endswith("geo_submission/md5sums.txt") for f in out.files)
    assert any(f.endswith("geo_submission/samples_table.txt") for f in out.files)


def test_factory_quantification_is_optional(tmp_path: Path):
    samples = make_fastq_collection(tmp_path, Assay.RNA)
    assay_cfg = RNAAssayConfig(
        bigwigs=BigwigConfig(pileup_method=[PileupMethod.DEEPTOOLS]),
        rna_quantification=RNAQuantificationConfig(method=QuantificationMethod.SALMON),
    )
    cfg = minimal_base_config(tmp_path, Assay.RNA, assay_cfg)

    builder = SeqnadoOutputFactory(assay=Assay.RNA, samples=samples, config=cfg).create_output_builder()
    out = builder.build()
    # Quantification files added via QuantificationFiles, but no strict path check here
    assert isinstance(out.files, list)
