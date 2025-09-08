from pathlib import Path
import pytest

from seqnado.config import GenomeConfig, STARIndex


def mk_star_genome(tmp_path: Path, name: str = "hg38") -> GenomeConfig:
    idx_dir = tmp_path / "star_index"
    idx_dir.mkdir()
    return GenomeConfig(name=name, index=STARIndex(prefix=idx_dir))


def test_genome_predicts_organism(tmp_path: Path):
    genome = mk_star_genome(tmp_path, name="mm10")
    assert genome.organism == "Mus musculus"


def test_genomeconfig_missing_fasta_raises_userfriendlyerror(tmp_path: Path):
    from seqnado.config import UserFriendlyError

    genome = mk_star_genome(tmp_path)
    with pytest.raises(UserFriendlyError) as ei:
        GenomeConfig(name=genome.name, index=genome.index, fasta=tmp_path / "no.fa")
    assert "Genome file not found" in str(ei.value)
