"""Tests for seqnado.inputs.bigwigs module."""

import pytest
from pathlib import Path

from seqnado import Assay
from seqnado.inputs.bigwigs import BigWigFile, BigWigCollection
from seqnado.inputs.core import Metadata


class TestBigWigFile:
    """Tests for BigWigFile class."""

    def test_create_bigwig_file(self, tmp_path):
        """Test creating a BigWigFile."""
        bw_path = tmp_path / "sample1.bigWig"
        bw_path.touch()
        
        bw = BigWigFile(path=bw_path)
        assert bw.path.is_absolute()
        assert bw.stem == "sample1"

    def test_bigwig_file_not_found_raises_error(self, tmp_path):
        """Test that missing bigWig file raises FileNotFoundError."""
        bw_path = tmp_path / "missing.bigWig"
        
        with pytest.raises(FileNotFoundError, match="bigWig file not found"):
            BigWigFile(path=bw_path)

    def test_infer_sample_id_rna_plus(self, tmp_path):
        """Test inferring sample ID from RNA plus strand file."""
        bw_path = tmp_path / "sample1_plus.bigWig"
        bw_path.touch()
        bw = BigWigFile(path=bw_path)
        
        sample_id = bw.infer_sample_id(Assay.RNA)
        assert sample_id == "sample1"

    def test_infer_sample_id_rna_minus(self, tmp_path):
        """Test inferring sample ID from RNA minus strand file."""
        bw_path = tmp_path / "sample1_minus.bigWig"
        bw_path.touch()
        bw = BigWigFile(path=bw_path)
        
        sample_id = bw.infer_sample_id(Assay.RNA)
        assert sample_id == "sample1"

    def test_infer_sample_id_rna_no_strand(self, tmp_path):
        """Test inferring sample ID from RNA file without strand suffix."""
        bw_path = tmp_path / "sample1.bigWig"
        bw_path.touch()
        bw = BigWigFile(path=bw_path)
        
        sample_id = bw.infer_sample_id(Assay.RNA)
        assert sample_id == "sample1"

    def test_infer_sample_id_chip(self, tmp_path):
        """Test inferring sample ID from ChIP file."""
        bw_path = tmp_path / "sample1.bigWig"
        bw_path.touch()
        bw = BigWigFile(path=bw_path)
        
        sample_id = bw.infer_sample_id(Assay.CHIP)
        assert sample_id == "sample1"

    def test_is_strand_specific_rna_plus(self, tmp_path):
        """Test strand-specificity detection for RNA plus strand."""
        bw_path = tmp_path / "sample1_plus.bigWig"
        bw_path.touch()
        bw = BigWigFile(path=bw_path)
        
        assert bw.is_strand_specific(Assay.RNA) is True

    def test_is_strand_specific_rna_minus(self, tmp_path):
        """Test strand-specificity detection for RNA minus strand."""
        bw_path = tmp_path / "sample1_minus.bigWig"
        bw_path.touch()
        bw = BigWigFile(path=bw_path)
        
        assert bw.is_strand_specific(Assay.RNA) is True

    def test_is_strand_specific_rna_no_strand(self, tmp_path):
        """Test strand-specificity detection for RNA without strand."""
        bw_path = tmp_path / "sample1.bigWig"
        bw_path.touch()
        bw = BigWigFile(path=bw_path)
        
        assert bw.is_strand_specific(Assay.RNA) is False

    def test_is_strand_specific_chip(self, tmp_path):
        """Test strand-specificity detection for ChIP."""
        bw_path = tmp_path / "sample1_plus.bigWig"
        bw_path.touch()
        bw = BigWigFile(path=bw_path)
        
        assert bw.is_strand_specific(Assay.CHIP) is False

    def test_strand_rna_plus(self, tmp_path):
        """Test getting strand for RNA plus file."""
        bw_path = tmp_path / "sample1_plus.bigWig"
        bw_path.touch()
        bw = BigWigFile(path=bw_path)
        
        strand = bw.strand(Assay.RNA)
        assert strand == "plus"

    def test_strand_rna_minus(self, tmp_path):
        """Test getting strand for RNA minus file."""
        bw_path = tmp_path / "sample1_minus.bigWig"
        bw_path.touch()
        bw = BigWigFile(path=bw_path)
        
        strand = bw.strand(Assay.RNA)
        assert strand == "minus"

    def test_strand_rna_no_strand(self, tmp_path):
        """Test getting strand for RNA file without strand suffix."""
        bw_path = tmp_path / "sample1.bigWig"
        bw_path.touch()
        bw = BigWigFile(path=bw_path)
        
        strand = bw.strand(Assay.RNA)
        assert strand is None

    def test_strand_chip(self, tmp_path):
        """Test getting strand for ChIP file returns None."""
        bw_path = tmp_path / "sample1.bigWig"
        bw_path.touch()
        bw = BigWigFile(path=bw_path)
        
        strand = bw.strand(Assay.CHIP)
        assert strand is None


class TestBigWigCollection:
    """Tests for BigWigCollection class."""

    def test_create_collection_chip(self, tmp_path):
        """Test creating a BigWigCollection for ChIP."""
        bw1 = tmp_path / "sample1.bigWig"
        bw2 = tmp_path / "sample2.bigWig"
        bw1.touch()
        bw2.touch()
        
        bigwig_files = [BigWigFile(path=bw1), BigWigFile(path=bw2)]
        metadata = [
            Metadata(sample_id="sample1", assay=Assay.CHIP),
            Metadata(sample_id="sample2", assay=Assay.CHIP),
        ]
        
        collection = BigWigCollection(
            assay=Assay.CHIP,
            bigwig_files=bigwig_files,
            metadata=metadata
        )
        
        assert collection.assay == Assay.CHIP
        assert len(collection.bigwig_files) == 2

    def test_primary_file_type(self, tmp_path):
        """Test primary_file_type property."""
        bw = tmp_path / "sample1.bigWig"
        bw.touch()
        
        collection = BigWigCollection(
            assay=Assay.CHIP,
            bigwig_files=[BigWigFile(path=bw)],
            metadata=[Metadata(sample_id="sample1", assay=Assay.CHIP)]
        )
        
        assert collection.primary_file_type == "bigwig"

    def test_sample_ids_chip(self, tmp_path):
        """Test sample_ids property for ChIP."""
        bw1 = tmp_path / "sample1.bigWig"
        bw2 = tmp_path / "sample2.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection(
            assay=Assay.CHIP,
            bigwig_files=[BigWigFile(path=bw1), BigWigFile(path=bw2)],
            metadata=[
                Metadata(sample_id="sample1", assay=Assay.CHIP),
                Metadata(sample_id="sample2", assay=Assay.CHIP),
            ]
        )
        
        assert collection.sample_ids == ["sample1", "sample2"]

    def test_sample_ids_rna_stranded(self, tmp_path):
        """Test sample_ids property for stranded RNA."""
        bw1 = tmp_path / "sample1_plus.bigWig"
        bw2 = tmp_path / "sample1_minus.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection(
            assay=Assay.RNA,
            bigwig_files=[BigWigFile(path=bw1), BigWigFile(path=bw2)],
            metadata=[
                Metadata(sample_id="sample1", assay=Assay.RNA),
                Metadata(sample_id="sample1", assay=Assay.RNA),
            ]
        )
        
        assert collection.sample_ids == ["sample1", "sample1"]

    def test_sample_names_chip(self, tmp_path):
        """Test sample_names property for ChIP (unique and sorted)."""
        bw1 = tmp_path / "sample1.bigWig"
        bw2 = tmp_path / "sample2.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection(
            assay=Assay.CHIP,
            bigwig_files=[BigWigFile(path=bw1), BigWigFile(path=bw2)],
            metadata=[
                Metadata(sample_id="sample1", assay=Assay.CHIP),
                Metadata(sample_id="sample2", assay=Assay.CHIP),
            ]
        )
        
        assert collection.sample_names == ["sample1", "sample2"]

    def test_sample_names_rna_stranded(self, tmp_path):
        """Test sample_names deduplicates stranded RNA samples."""
        bw1 = tmp_path / "sample1_plus.bigWig"
        bw2 = tmp_path / "sample1_minus.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection(
            assay=Assay.RNA,
            bigwig_files=[BigWigFile(path=bw1), BigWigFile(path=bw2)],
            metadata=[
                Metadata(sample_id="sample1", assay=Assay.RNA),
                Metadata(sample_id="sample1", assay=Assay.RNA),
            ]
        )
        
        assert collection.sample_names == ["sample1"]

    def test_tracks_chip(self, tmp_path):
        """Test tracks property for ChIP samples."""
        bw1 = tmp_path / "sample1.bigWig"
        bw2 = tmp_path / "sample2.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection(
            assay=Assay.CHIP,
            bigwig_files=[BigWigFile(path=bw1), BigWigFile(path=bw2)],
            metadata=[
                Metadata(sample_id="sample1", assay=Assay.CHIP),
                Metadata(sample_id="sample2", assay=Assay.CHIP),
            ]
        )
        
        tracks = collection.tracks
        assert "sample1" in tracks
        assert "sample2" in tracks
        assert tracks["sample1"] == [bw1]
        assert tracks["sample2"] == [bw2]

    def test_tracks_rna_stranded(self, tmp_path):
        """Test tracks property groups stranded RNA files."""
        bw1 = tmp_path / "sample1_plus.bigWig"
        bw2 = tmp_path / "sample1_minus.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection(
            assay=Assay.RNA,
            bigwig_files=[BigWigFile(path=bw1), BigWigFile(path=bw2)],
            metadata=[
                Metadata(sample_id="sample1", assay=Assay.RNA),
                Metadata(sample_id="sample1", assay=Assay.RNA),
            ]
        )
        
        tracks = collection.tracks
        assert "sample1" in tracks
        assert len(tracks["sample1"]) == 2
        # Check that plus comes before minus
        assert "plus" in tracks["sample1"][0].stem
        assert "minus" in tracks["sample1"][1].stem

    def test_get_file_paths_bigwig(self, tmp_path):
        """Test get_file_paths with bigwig kind."""
        bw1 = tmp_path / "sample1.bigWig"
        bw2 = tmp_path / "sample2.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection(
            assay=Assay.CHIP,
            bigwig_files=[BigWigFile(path=bw1), BigWigFile(path=bw2)],
            metadata=[
                Metadata(sample_id="sample1", assay=Assay.CHIP),
                Metadata(sample_id="sample2", assay=Assay.CHIP),
            ]
        )
        
        paths = collection.get_file_paths("bigwig")
        assert len(paths) == 2
        assert bw1 in paths
        assert bw2 in paths

    def test_get_file_paths_none(self, tmp_path):
        """Test get_file_paths with None returns all bigwigs."""
        bw1 = tmp_path / "sample1.bigWig"
        bw1.touch()
        
        collection = BigWigCollection(
            assay=Assay.CHIP,
            bigwig_files=[BigWigFile(path=bw1)],
            metadata=[Metadata(sample_id="sample1", assay=Assay.CHIP)]
        )
        
        paths = collection.get_file_paths(None)
        assert paths == [bw1]

    def test_get_file_paths_invalid_kind(self, tmp_path):
        """Test get_file_paths raises error for invalid kind."""
        bw1 = tmp_path / "sample1.bigWig"
        bw1.touch()
        
        collection = BigWigCollection(
            assay=Assay.CHIP,
            bigwig_files=[BigWigFile(path=bw1)],
            metadata=[Metadata(sample_id="sample1", assay=Assay.CHIP)]
        )
        
        with pytest.raises(ValueError, match="Unsupported file kind 'bam'"):
            collection.get_file_paths("bam")

    def test_from_files_chip(self, tmp_path):
        """Test creating BigWigCollection from files list for ChIP."""
        bw1 = tmp_path / "sample1.bigWig"
        bw2 = tmp_path / "sample2.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection.from_files(
            assay=Assay.CHIP,
            files=[bw1, bw2]
        )
        
        assert len(collection.bigwig_files) == 2
        assert len(collection.metadata) == 2
        assert collection.sample_names == ["sample1", "sample2"]

    def test_from_files_rna_stranded(self, tmp_path):
        """Test creating BigWigCollection from stranded RNA files."""
        bw1 = tmp_path / "sample1_plus.bigWig"
        bw2 = tmp_path / "sample1_minus.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection.from_files(
            assay=Assay.RNA,
            files=[bw1, bw2]
        )
        
        assert len(collection.bigwig_files) == 2
        assert collection.sample_names == ["sample1"]

    def test_from_files_with_metadata_function(self, tmp_path):
        """Test from_files with metadata function."""
        bw1 = tmp_path / "sample1.bigWig"
        bw1.touch()
        
        def metadata_fn(sample_id: str) -> Metadata:
            return Metadata(assay=Assay.CHIP, scaling_group="batch1")
        
        collection = BigWigCollection.from_files(
            assay=Assay.CHIP,
            files=[bw1],
            metadata=metadata_fn
        )
        
        assert collection.metadata[0].scaling_group == "batch1"

    def test_from_files_with_metadata_kwargs(self, tmp_path):
        """Test from_files with metadata kwargs."""
        bw1 = tmp_path / "sample1.bigWig"
        bw1.touch()
        
        collection = BigWigCollection.from_files(
            assay=Assay.CHIP,
            files=[bw1],
            scaling_group="batch2"
        )
        
        assert collection.metadata[0].scaling_group == "batch2"

    def test_from_files_empty_raises_error(self):
        """Test from_files with no files raises ValueError."""
        with pytest.raises(ValueError, match="No bigWig files provided"):
            BigWigCollection.from_files(assay=Assay.CHIP, files=[])

    def test_from_directory_bigwig_extension(self, tmp_path):
        """Test from_directory finds .bigWig files."""
        bw1 = tmp_path / "sample1.bigWig"
        bw2 = tmp_path / "sample2.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection.from_directory(
            assay=Assay.CHIP,
            directory=tmp_path
        )
        
        assert len(collection.bigwig_files) == 2

    def test_from_directory_bw_extension(self, tmp_path):
        """Test from_directory finds .bw files."""
        bw1 = tmp_path / "sample1.bw"
        bw2 = tmp_path / "sample2.bw"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection.from_directory(
            assay=Assay.CHIP,
            directory=tmp_path
        )
        
        assert len(collection.bigwig_files) == 2

    def test_from_directory_recursive(self, tmp_path):
        """Test from_directory finds files in subdirectories."""
        subdir = tmp_path / "subdir"
        subdir.mkdir()
        bw1 = subdir / "sample1.bigWig"
        bw1.touch()
        
        collection = BigWigCollection.from_directory(
            assay=Assay.CHIP,
            directory=tmp_path
        )
        
        assert len(collection.bigwig_files) == 1

    def test_from_directory_no_files_raises_error(self, tmp_path):
        """Test from_directory raises error when no files found."""
        with pytest.raises(FileNotFoundError, match="No bigWig files found"):
            BigWigCollection.from_directory(
                assay=Assay.CHIP,
                directory=tmp_path
            )

    def test_from_directory_custom_glob_patterns(self, tmp_path):
        """Test from_directory with custom glob patterns."""
        bw1 = tmp_path / "sample1.custom"
        bw1.touch()
        
        collection = BigWigCollection.from_directory(
            assay=Assay.CHIP,
            directory=tmp_path,
            glob_patterns=("*.custom",)
        )
        
        assert len(collection.bigwig_files) == 1

    def test_to_dataframe_chip(self, tmp_path):
        """Test to_dataframe for ChIP samples."""
        bw1 = tmp_path / "sample1.bigWig"
        bw2 = tmp_path / "sample2.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection.from_files(
            assay=Assay.CHIP,
            files=[bw1, bw2],
            scaling_group="batch1"
        )
        
        df = collection.to_dataframe()
        
        assert len(df) == 2
        assert "sample_id" in df.columns
        assert "bigwig" in df.columns
        assert df["scaling_group"].iloc[0] == "batch1"

    def test_to_dataframe_rna_stranded(self, tmp_path):
        """Test to_dataframe for stranded RNA samples."""
        bw1 = tmp_path / "sample1_plus.bigWig"
        bw2 = tmp_path / "sample1_minus.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection.from_files(
            assay=Assay.RNA,
            files=[bw1, bw2]
        )
        
        df = collection.to_dataframe()
        
        assert len(df) == 1
        assert "sample_id" in df.columns
        assert "bigwig_plus" in df.columns
        assert "bigwig_minus" in df.columns
        assert df["bigwig_plus"].iloc[0] == bw1
        assert df["bigwig_minus"].iloc[0] == bw2

    def test_to_dataframe_sorted_by_sample_id(self, tmp_path):
        """Test to_dataframe returns sorted DataFrame."""
        bw1 = tmp_path / "sample2.bigWig"
        bw2 = tmp_path / "sample1.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection.from_files(
            assay=Assay.CHIP,
            files=[bw1, bw2]
        )
        
        df = collection.to_dataframe()
        
        # Check sorted order by looking at values
        sample_ids = df["sample_id"].tolist()
        assert sample_ids == ["sample1", "sample2"]

    def test_symlink_bigwig_files(self, tmp_path):
        """Test symlinking bigWig files to output directory."""
        source_dir = tmp_path / "source"
        source_dir.mkdir()
        bw1 = source_dir / "sample1.bigWig"
        bw1.touch()
        
        collection = BigWigCollection.from_files(
            assay=Assay.CHIP,
            files=[bw1]
        )
        
        output_dir = tmp_path / "output"
        collection.symlink_bigwig_files(output_dir)
        
        link = output_dir / "sample1.bigWig"
        assert link.exists()
        assert link.is_symlink()
        assert link.resolve() == bw1.resolve()

    def test_symlink_bigwig_files_creates_directory(self, tmp_path):
        """Test that symlink_bigwig_files creates output directory."""
        bw1 = tmp_path / "sample1.bigWig"
        bw1.touch()
        
        collection = BigWigCollection.from_files(
            assay=Assay.CHIP,
            files=[bw1]
        )
        
        output_dir = tmp_path / "new_output"
        collection.symlink_bigwig_files(output_dir)
        
        assert output_dir.exists()
        assert (output_dir / "sample1.bigWig").exists()

    def test_symlink_bigwig_files_skips_existing(self, tmp_path):
        """Test that symlink_bigwig_files doesn't overwrite existing links."""
        bw1 = tmp_path / "sample1.bigWig"
        bw1.touch()
        
        collection = BigWigCollection.from_files(
            assay=Assay.CHIP,
            files=[bw1]
        )
        
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        link = output_dir / "sample1.bigWig"
        link.symlink_to(bw1)
        
        # Should not raise error
        collection.symlink_bigwig_files(output_dir)
        assert link.exists()
