from pathlib import Path

from seqnado import Assay, PileupMethod
from seqnado.config.configs import BigwigConfig, GenomeConfig, STARIndex
from seqnado.config.core import ATACAssayConfig, SeqnadoConfig
from seqnado.inputs.fastq import FastqCollection, FastqFile, FastqSet
from seqnado.inputs.grouping import SampleGroup, SampleGroupings, SampleGroups
from seqnado.outputs.core import SeqnadoOutputBuilder, SeqnadoOutputFactory


def _minimal_config(tmp: Path) -> SeqnadoConfig:
    star = tmp / "star"
    star.mkdir()
    genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
    assay_cfg = ATACAssayConfig(
        bigwigs=BigwigConfig(pileup_method=[PileupMethod.DEEPTOOLS])
    )
    return SeqnadoConfig(
        assay=Assay.ATAC,
        project=dict(name="p"),
        genome=genome,
        metadata=tmp / "m.csv",
        assay_config=assay_cfg,
    )


def _small_collection(tmp: Path) -> FastqCollection:
    r1_path = tmp / "s1_R1.fastq.gz"
    r1_path.write_text("@r\nN\n+\n#\n")
    r2_path = tmp / "s1_R2.fastq.gz"
    r2_path.write_text("@r\nN\n+\n#\n")
    r1 = FastqFile(path=r1_path)
    r2 = FastqFile(path=r2_path)
    fs = FastqSet(sample_id="s1", r1=r1, r2=r2)
    return FastqCollection(assay=Assay.ATAC, metadata=[], fastq_sets=[fs])


def test_output_builder_bigwigs_only(tmp_path: Path):
    cfg = _minimal_config(tmp_path)
    samples = _small_collection(tmp_path)

    builder = SeqnadoOutputBuilder(Assay.ATAC, samples, cfg)
    builder.add_qc_files()
    builder.add_individual_bigwig_files()
    out = builder.build()

    assert any("seqnado_report.html" in f for f in out.files)
    assert any(f.endswith(".bigWig") for f in out.files)


def test_output_factory_with_consensus_groups(tmp_path: Path):
    cfg = _minimal_config(tmp_path)
    samples = _small_collection(tmp_path)
    groups = SampleGroupings(
        groupings={
            "consensus": SampleGroups(groups=[SampleGroup(name="cons", samples=["s1"])])
        }
    )
    factory = SeqnadoOutputFactory(Assay.ATAC, samples, cfg, sample_groupings=groups)
    out = factory.create_output_builder().build()
    assert out.files, "Expected some files listed from factory"
