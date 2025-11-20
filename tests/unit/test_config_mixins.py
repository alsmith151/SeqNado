from pathlib import Path

from seqnado import Assay
from seqnado.config.configs import (
    GenomeConfig, STARIndex, BigwigConfig, PlottingConfig, UCSCHubConfig,
)
from seqnado.config.core import SeqnadoConfig, ATACAssayConfig


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