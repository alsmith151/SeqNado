from pathlib import Path

from seqnado import Assay
from seqnado.config.configs import (
    BigwigConfig,
    GenomeConfig,
    PlottingConfig,
    STARIndex,
    UCSCHubConfig,
)
from seqnado.config.core import ATACAssayConfig, SeqnadoConfig


def test_common_computed_flags(tmp_path: Path):
    # Minimal STAR index dir
    star_dir = tmp_path / "star"
    star_dir.mkdir()

    # Create parent directory for hub (validator checks parent exists)
    hub_parent = tmp_path
    hub_parent.mkdir(exist_ok=True)

    genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star_dir))
    assay_cfg = ATACAssayConfig(
        bigwigs=BigwigConfig(pileup_method=[]),
        plotting=PlottingConfig(coordinates=None),
        ucsc_hub=UCSCHubConfig(directory=str(tmp_path / "hub")),
    )

    cfg = SeqnadoConfig(
        assay=Assay.ATAC,
        project=dict(name="p"),
        genome=genome,
        metadata=tmp_path / "m.csv",
        assay_config=assay_cfg,
    )

    assert cfg.assay_config.create_bigwigs is True
    assert cfg.assay_config.create_ucsc_hub is True
    assert cfg.assay_config.plot_with_plotnado is True


class TestSeqnadoConfigComputedFields:
    """Tests for SeqnadoConfig computed fields and validators."""

    def test_organism_computed_field(self, tmp_path):
        """Test organism computed field."""
        from seqnado import Organism
        
        star_dir = tmp_path / "star"
        star_dir.mkdir()
        
        genome = GenomeConfig(
            name="hg38",
            index=STARIndex(prefix=star_dir),
            organism=Organism.HUMAN.value,
        )
        
        cfg = SeqnadoConfig(
            assay=Assay.RNA,
            project=dict(name="test"),
            genome=genome,
            metadata=tmp_path / "m.csv",
        )
        
        assert cfg.organism == Organism.HUMAN.value

    def test_shift_for_tn5_insertion_true(self, tmp_path):
        """Test shift_for_tn5_insertion when tn5_shift is True."""
        star_dir = tmp_path / "star"
        star_dir.mkdir()
        
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star_dir))
        assay_cfg = ATACAssayConfig(tn5_shift=True)
        
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project=dict(name="test"),
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )
        
        assert cfg.shift_for_tn5_insertion is True

    def test_shift_for_tn5_insertion_false(self, tmp_path):
        """Test shift_for_tn5_insertion when tn5_shift is False."""
        star_dir = tmp_path / "star"
        star_dir.mkdir()
        
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star_dir))
        assay_cfg = ATACAssayConfig(tn5_shift=False)
        
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project=dict(name="test"),
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )
        
        assert cfg.shift_for_tn5_insertion is False

    def test_mcc_viewpoints_when_absent(self, tmp_path):
        """Test mcc_viewpoints returns empty string when not MCC."""
        star_dir = tmp_path / "star"
        star_dir.mkdir()
        
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star_dir))
        
        cfg = SeqnadoConfig(
            assay=Assay.RNA,
            project=dict(name="test"),
            genome=genome,
            metadata=tmp_path / "m.csv",
        )
        
        assert cfg.mcc_viewpoints == ""
