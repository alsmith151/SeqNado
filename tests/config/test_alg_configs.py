from pathlib import Path
import pytest

from seqnado import Assay, SNPCallingMethod, SpikeInMethod
from seqnado.config import (
    SeqnadoConfig,
    ProjectConfig,
    GenomeConfig,
    STARIndex,
    ATACAssayConfig,
    ChIPAssayConfig,
    MCCAssayConfig,
    PeakCallingConfig,
    MLDatasetConfig,
    MCCConfig,
    SNPCallingConfig,
)


def mk_star_genome(tmp_path: Path, name: str = "hg38") -> GenomeConfig:
    d = tmp_path / "star"
    d.mkdir()
    return GenomeConfig(name=name, index=STARIndex(prefix=d))


def test_snp_calling_annotate_flag():
    s = SNPCallingConfig(method=SNPCallingMethod.BCFTOOLS, snp_database="dbsnp")
    assert s.annotate_snps is True


def test_snp_calling_annotate_flag_false_when_no_db():
    s = SNPCallingConfig(method=SNPCallingMethod.BCFTOOLS, snp_database=None)
    assert s.annotate_snps is False


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


def test_mcc_config_invalid_resolutions_raises():
    with pytest.raises(ValueError):
        MCCConfig(viewpoints=Path("/non/existent"), resolutions=[0])


def test_atac_assay_flags_and_seqnado_shift(tmp_path: Path):
    genome = mk_star_genome(tmp_path)
    atac = ATACAssayConfig(
        tn5_shift=True,
        bigwigs=None,
        plotting=None,
        dataset_for_ml=MLDatasetConfig(binsize=10),
        peak_calling=PeakCallingConfig(method=[]),
    )
    # Computed flags
    assert atac.create_dataset
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
    from seqnado.config import SpikeInConfig

    cfg = ChIPAssayConfig(spikein=SpikeInConfig(method=SpikeInMethod.ORLANDO), peak_calling=PeakCallingConfig(method=[]))
    assert getattr(cfg, "has_spikein") is True
    assert getattr(cfg, "call_peaks") is True

