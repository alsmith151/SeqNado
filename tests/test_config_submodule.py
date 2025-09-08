from pathlib import Path
import pytest

from seqnado import Assay, GenomicCoordinate, SNPCallingMethod, SpikeInMethod
from seqnado.config import (
    SeqnadoConfig,
    ATACAssayConfig,
    ChIPAssayConfig,
    CATAssayConfig,
    RNAAssayConfig,
    SNPAssayConfig,
    MCCAssayConfig,
    MethylationAssayConfig,
    STARIndex,
    GenomeConfig,
    BigwigConfig,
    ProjectConfig,
    PCRDuplicatesConfig,
    PlottingConfig,
    PeakCallingConfig,
    SpikeInConfig,
    UCSCHubConfig,
    RNAQuantificationConfig,
    SNPCallingConfig,
    MCCConfig,
    MethylationConfig,
    MLDatasetConfig,
    UserFriendlyError,
)


def mk_star_genome(tmp_path: Path, name: str = "hg38") -> GenomeConfig:
    idx_dir = tmp_path / "star_index"
    idx_dir.mkdir()
    return GenomeConfig(name=name, index=STARIndex(prefix=idx_dir))


def test_genome_predicts_organism(tmp_path: Path):
    genome = mk_star_genome(tmp_path, name="mm10")
    assert genome.organism == "Mus musculus"


def test_star_index_requires_existing_dir(tmp_path: Path):
    with pytest.raises(ValueError):
        STARIndex(prefix=tmp_path / "missing")
    # Valid when dir exists
    d = tmp_path / "ok"
    d.mkdir()
    s = STARIndex(prefix=d)
    assert isinstance(s, STARIndex)
    assert s.type == "STAR"


def test_genomic_coordinate_parsing_and_validation():
    gc = GenomicCoordinate.from_string("chr2:10-20")
    assert gc.chromosome == "chr2" and gc.start == 10 and gc.end == 20
    with pytest.raises(ValueError):
        GenomicCoordinate(chromosome="chr1", start=10, end=5)


def test_ucsc_hub_validation_and_fields(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    # Ensure parent dir exists for default path used in factory
    monkeypatch.chdir(tmp_path)
    (tmp_path / "seqnado_output").mkdir()

    # Construct explicitly to avoid default factory bug; check fields validate
    hub = UCSCHubConfig(
        directory=str(tmp_path / "hub"),
        name="seqnado_hub",
        genome="hg38",
        email="a@b.com",
        default_position=GenomicCoordinate(chromosome="chr1", start=100, end=200),
        subgroup_by=["method", "norm", "strand"],
        overlay_by=["samplename", "method", "norm"],
    )
    assert "strand" in hub.subgroup_by
    assert set(["samplename", "method", "norm"]).issubset(hub.overlay_by)

    # Name validation and directory parent existence
    with pytest.raises(ValueError):
        UCSCHubConfig(
            directory=str(tmp_path / "missing" / "hub"),
            name="ok",
            default_position=GenomicCoordinate.from_string("chr1:100-200"),
        )
    ok = UCSCHubConfig(directory=str(tmp_path / "hub"), name="seqnado_hub")
    assert ok.name == "seqnado_hub"
    with pytest.raises(ValueError):
        UCSCHubConfig(directory=str(tmp_path / "hub"), name="bad name")


def test_mldataset_requires_one_of_binsize_or_regions(tmp_path: Path):
    with pytest.raises(ValueError):
        MLDatasetConfig()
    # binsize must be positive
    with pytest.raises(ValueError):
        MLDatasetConfig(binsize=-1)
    # regions_bed path must exist when provided
    with pytest.raises(ValueError):
        MLDatasetConfig(binsize=100, regions_bed=tmp_path / "nope.bed")
    # happy path
    ok = MLDatasetConfig(binsize=50)
    assert ok.binsize == 50


def test_snp_calling_annotate_flag():
    s = SNPCallingConfig(method=SNPCallingMethod.BCFTOOLS, snp_database="dbsnp")
    assert s.annotate_snps is True


def test_mcc_config_validates_and_seqnado_computed_viewpoints(tmp_path: Path):
    vp = tmp_path / "vp.bed"
    vp.write_text("chr1\t100\t200\tvp1\n")
    mcc = MCCConfig(viewpoints=vp)
    assert mcc.resolutions == [100]

    genome = mk_star_genome(tmp_path)
    cfg = SeqnadoConfig(
        assay=Assay.MCC,
        project=ProjectConfig(name="proj"),
        genome=genome,
        metadata=tmp_path / "meta.csv",
        assay_config=MCCAssayConfig(mcc=mcc),
    )
    assert cfg.mcc_viewpoints == str(vp)


def test_atac_assay_flags_and_seqnado_shift(tmp_path: Path):
    genome = mk_star_genome(tmp_path)
    atac = ATACAssayConfig(
        tn5_shift=True,
        bigwigs=BigwigConfig(),
        plotting=PlottingConfig(),
        dataset_for_ml=MLDatasetConfig(binsize=10),
        peak_calling=PeakCallingConfig(),
    )
    # Computed flags
    assert atac.create_bigwigs and atac.plot_with_plotnado and atac.create_dataset
    assert atac.call_peaks is True

    seq = SeqnadoConfig(
        assay=Assay.ATAC,
        project=ProjectConfig(name="p"),
        genome=genome,
        metadata=tmp_path / "meta.csv",
        assay_config=atac,
    )
    assert seq.shift_for_tn5_insertion is True


def test_chip_has_spikein_and_calls_peaks():
    cfg = ChIPAssayConfig(spikein=SpikeInConfig(method=SpikeInMethod.ORLANDO), peak_calling=PeakCallingConfig())
    assert getattr(cfg, "has_spikein") is True
    assert getattr(cfg, "call_peaks") is True


def test_qc_config_defaults_true():
    from seqnado.config.configs import QCConfig

    qc = QCConfig()
    assert qc.run_fastq_screen and qc.calculate_library_size and qc.calculate_fraction_of_reads_in_peaks


# -----------------------
# More complex, example-driven tests
# -----------------------

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


def test_genomeconfig_missing_fasta_raises_userfriendlyerror(tmp_path: Path):
    genome = mk_star_genome(tmp_path)
    with pytest.raises(UserFriendlyError) as ei:
        GenomeConfig(name=genome.name, index=genome.index, fasta=tmp_path / "no.fa")
    assert "Genome file not found" in str(ei.value)


def test_seqnadoconfig_assay_config_dict_autocoerce(tmp_path: Path):
    genome = mk_star_genome(tmp_path)
    cfg = SeqnadoConfig(
        assay=Assay.ATAC,
        project=ProjectConfig(name="proj"),
        genome=genome,
        metadata=tmp_path / "meta.csv",
        assay_config={"tn5_shift": True},
    )
    assert isinstance(cfg.assay_config, ATACAssayConfig)
    assert cfg.shift_for_tn5_insertion is True


def test_snp_calling_annotate_flag_false_when_no_db():
    s = SNPCallingConfig(method=SNPCallingMethod.BCFTOOLS, snp_database=None)
    assert s.annotate_snps is False


def test_mcc_config_invalid_resolutions_raises():
    with pytest.raises(ValueError):
        MCCConfig(viewpoints=Path("/non/existent"), resolutions=[0])


def test_bigwig_config_enum_parsing_and_binsize():
    from seqnado import PileupMethod

    bw = BigwigConfig(pileup_method=["deeptools", PileupMethod.HOMER], binsize=25)
    assert bw.pileup_method and set(bw.pileup_method) == {PileupMethod.DEEPTOOLS, PileupMethod.HOMER}
    assert bw.binsize == 25


def test_rna_quantification_salmon_index_path_validator(tmp_path: Path):
    with pytest.raises(ValueError):
        RNAQuantificationConfig(method=None, salmon_index=str(tmp_path / "missing"))  # type: ignore[arg-type]