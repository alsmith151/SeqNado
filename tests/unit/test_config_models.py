from pathlib import Path

import pytest

from seqnado.config.configs import (
    BowtieIndex,
    GenomeConfig,
    STARIndex,
    UserFriendlyError,
)


def test_bowtie_index_files_validation(tmp_path: Path):
    # Create fake bowtie2 index files with typical suffixes
    idx_prefix = tmp_path / "idx/genome"
    idx_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in [".1.bt2", ".2.bt2", ".rev.1.bt2"]:
        (Path(str(idx_prefix) + suffix)).write_text("")

    cfg = GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(idx_prefix)))
    assert isinstance(cfg.index, BowtieIndex)
    assert cfg.index.files, "Expected to find bowtie index files"


def test_star_index_validation(tmp_path: Path):
    idx_dir = tmp_path / "star"
    idx_dir.mkdir(parents=True)
    (idx_dir / "SA").write_text("")  # dummy content
    cfg = GenomeConfig(name="hg38", index=STARIndex(prefix=idx_dir))
    assert isinstance(cfg.index, STARIndex)


def test_missing_optional_paths_raise_userfriendlyerror(tmp_path: Path):
    # fasta points to a missing file -> should raise our friendly error
    with pytest.raises(UserFriendlyError):
        GenomeConfig(
            name="hg38", index=STARIndex(prefix=tmp_path), fasta=tmp_path / "missing.fa"
        )


class TestGenomeConfigValidation:
    """Tests for GenomeConfig path validation and computed fields."""

    def test_genome_config_with_all_paths(self, tmp_path):
        """Test GenomeConfig with all optional paths provided."""
        idx_dir = tmp_path / "star"
        idx_dir.mkdir()
        (idx_dir / "SA").write_text("")
        
        fasta = tmp_path / "genome.fa"
        fasta.touch()
        gtf = tmp_path / "genes.gtf"
        gtf.touch()
        blacklist = tmp_path / "blacklist.bed"
        blacklist.touch()
        
        config = GenomeConfig(
            name="hg38",
            index=STARIndex(prefix=idx_dir),
            fasta=fasta,
            gtf=gtf,
            blacklist=blacklist,
        )
        
        assert config.name == "hg38"
        assert config.fasta == fasta
        assert config.gtf == gtf
        assert config.blacklist == blacklist

    def test_genome_config_missing_gtf(self, tmp_path):
        """Test GenomeConfig raises error for missing GTF."""
        idx_dir = tmp_path / "star"
        idx_dir.mkdir()
        (idx_dir / "SA").write_text("")
        
        with pytest.raises(UserFriendlyError):
            GenomeConfig(
                name="hg38",
                index=STARIndex(prefix=idx_dir),
                gtf=tmp_path / "missing.gtf",
            )

    def test_genome_config_missing_blacklist(self, tmp_path):
        """Test GenomeConfig raises error for missing blacklist."""
        idx_dir = tmp_path / "star"
        idx_dir.mkdir()
        (idx_dir / "SA").write_text("")
        
        with pytest.raises(UserFriendlyError):
            GenomeConfig(
                name="hg38",
                index=STARIndex(prefix=idx_dir),
                blacklist=tmp_path / "missing.bed",
            )

    def test_bowtie_index_prefix_normalization(self, tmp_path):
        """Test BowtieIndex normalizes prefix paths."""
        idx_prefix = tmp_path / "idx/genome"
        idx_prefix.parent.mkdir(parents=True, exist_ok=True)
        for suffix in [".1.bt2", ".2.bt2", ".rev.1.bt2"]:
            (Path(str(idx_prefix) + suffix)).write_text("")
        
        index = BowtieIndex(prefix=str(idx_prefix))
        assert index.prefix == str(idx_prefix)
        assert len(index.files) > 0

    def test_star_index_directory_validation(self, tmp_path):
        """Test STARIndex validates directory exists."""
        idx_dir = tmp_path / "star"
        idx_dir.mkdir()
        
        index = STARIndex(prefix=idx_dir)
        assert index.prefix == idx_dir

    def test_star_index_missing_directory_fails(self, tmp_path):
        """Test STARIndex fails for missing directory."""
        with pytest.raises(ValueError, match="does not exist"):
            STARIndex(prefix=tmp_path / "missing")

    def test_bowtie_index_missing_directory_fails(self, tmp_path):
        """Test BowtieIndex fails for missing directory."""
        with pytest.raises(ValueError, match="does not exist"):
            BowtieIndex(prefix=str(tmp_path / "missing/index"))


