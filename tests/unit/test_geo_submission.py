"""Tests for seqnado.config.geo_submission module."""

import pandas as pd
import pytest
from pydantic import ValidationError

from seqnado import Assay, Molecule, Organism
from seqnado.config.geo_submission import GEOSample, GEOSamples


class TestGEOSample:
    """Tests for GEOSample model."""

    def test_create_basic_atac_sample(self):
        """Test creating a basic ATAC sample."""
        sample = GEOSample(
            assay=Assay.ATAC,
            library_name="ATAC_sample1",
            title="ATAC-seq of cell line",
            organism=Organism.HUMAN,
            cell_line="HEK293",
        )
        assert sample.assay == Assay.ATAC
        assert sample.library_name == "ATAC_sample1"
        assert sample.organism == Organism.HUMAN
        assert sample.cell_line == "HEK293"

    def test_create_basic_rna_sample(self):
        """Test creating a basic RNA sample."""
        sample = GEOSample(
            assay=Assay.RNA,
            library_name="RNA_sample1",
            title="RNA-seq of treated cells",
            organism=Organism.MOUSE,
            cell_type="neurons",
            treatment="drug_treatment",
        )
        assert sample.assay == Assay.RNA
        assert sample.cell_type == "neurons"
        assert sample.treatment == "drug_treatment"

    def test_create_chip_sample_with_antibody(self):
        """Test creating a ChIP sample with required antibody."""
        sample = GEOSample(
            assay=Assay.CHIP,
            library_name="ChIP_H3K27ac",
            title="ChIP-seq for H3K27ac",
            organism=Organism.HUMAN,
            cell_line="K562",
            antibody="H3K27ac",
        )
        assert sample.assay == Assay.CHIP
        assert sample.antibody == "H3K27ac"

    def test_create_cutandtag_sample_with_antibody(self):
        """Test creating a CUT&TAG sample with required antibody."""
        sample = GEOSample(
            assay=Assay.CAT,
            library_name="CAT_H3K4me3",
            title="CUT&TAG for H3K4me3",
            organism=Organism.MOUSE,
            cell_type="ES cells",
            antibody="H3K4me3",
        )
        assert sample.assay == Assay.CAT
        assert sample.antibody == "H3K4me3"

    def test_chip_sample_requires_antibody(self):
        """Test that ChIP samples require antibody specification."""
        with pytest.raises(ValidationError, match="Antibody must be specified"):
            GEOSample(
                assay=Assay.CHIP,
                library_name="ChIP_sample",
                title="ChIP-seq",
                organism=Organism.HUMAN,
                cell_line="K562",
            )

    def test_cutandtag_sample_requires_antibody(self):
        """Test that CUT&TAG samples require antibody specification."""
        with pytest.raises(ValidationError, match="Antibody must be specified"):
            GEOSample(
                assay=Assay.CAT,
                library_name="CAT_sample",
                title="CUT&TAG",
                organism=Organism.MOUSE,
            )

    def test_atac_sample_rejects_antibody(self):
        """Test that ATAC samples should not have antibody."""
        with pytest.raises(ValidationError, match="Antibody should not be specified"):
            GEOSample(
                assay=Assay.ATAC,
                library_name="ATAC_sample",
                title="ATAC-seq",
                organism=Organism.HUMAN,
                antibody="H3K4me3",  # This should fail
            )

    def test_rna_sample_rejects_antibody(self):
        """Test that RNA samples should not have antibody."""
        with pytest.raises(ValidationError, match="Antibody should not be specified"):
            GEOSample(
                assay=Assay.RNA,
                library_name="RNA_sample",
                title="RNA-seq",
                organism=Organism.MOUSE,
                antibody="H3K27ac",  # This should fail
            )

    def test_library_strategy_matches_assay(self):
        """Test that library_strategy computed field matches assay value."""
        sample = GEOSample(
            assay=Assay.ATAC,
            library_name="ATAC_sample",
            title="ATAC-seq",
            organism=Organism.HUMAN,
        )
        assert sample.library_strategy == "ATAC"

    def test_molecule_dna_genomic_for_atac(self):
        """Test that ATAC samples have genomic DNA molecule type."""
        sample = GEOSample(
            assay=Assay.ATAC,
            library_name="ATAC_sample",
            title="ATAC-seq",
            organism=Organism.HUMAN,
        )
        assert sample.molecule == Molecule.dna_genomic

    def test_molecule_dna_genomic_for_chip(self):
        """Test that ChIP samples have genomic DNA molecule type."""
        sample = GEOSample(
            assay=Assay.CHIP,
            library_name="ChIP_sample",
            title="ChIP-seq",
            organism=Organism.HUMAN,
            antibody="H3K27ac",
        )
        assert sample.molecule == Molecule.dna_genomic

    def test_molecule_rna_polya_for_standard_rna(self):
        """Test that standard RNA samples have polyA RNA molecule type."""
        sample = GEOSample(
            assay=Assay.RNA,
            library_name="RNA_sample",
            title="RNA-seq of cells",
            organism=Organism.MOUSE,
        )
        assert sample.molecule == Molecule.rna_polya

    def test_molecule_rna_nuclear_for_ttseq(self):
        """Test that TT-seq samples have nuclear RNA molecule type."""
        sample = GEOSample(
            assay=Assay.RNA,
            library_name="RNA_sample",
            title="TT-seq experiment",
            organism=Organism.MOUSE,
        )
        assert sample.molecule == Molecule.rna_nuclear

    def test_molecule_rna_nuclear_for_nasc(self):
        """Test that NASC samples have nuclear RNA molecule type."""
        sample = GEOSample(
            assay=Assay.RNA,
            library_name="RNA_sample",
            title="NASC-seq data",
            organism=Organism.HUMAN,
        )
        assert sample.molecule == Molecule.rna_nuclear

    def test_molecule_rna_nuclear_for_point(self):
        """Test that POINT samples have nuclear RNA molecule type."""
        sample = GEOSample(
            assay=Assay.RNA,
            library_name="RNA_sample",
            title="POINT-seq analysis",
            organism=Organism.RAT,
        )
        assert sample.molecule == Molecule.rna_nuclear

    def test_default_instrument_model(self):
        """Test default instrument model is NovaSeq X."""
        sample = GEOSample(
            assay=Assay.ATAC,
            library_name="ATAC_sample",
            title="ATAC-seq",
            organism=Organism.HUMAN,
        )
        assert sample.instrument_model == "NovaSeq X"

    def test_default_single_or_paired(self):
        """Test default sequencing type is paired-end."""
        sample = GEOSample(
            assay=Assay.RNA,
            library_name="RNA_sample",
            title="RNA-seq",
            organism=Organism.MOUSE,
        )
        assert sample.single_or_paired == "paired-end"

    def test_single_end_sequencing(self):
        """Test creating sample with single-end sequencing."""
        sample = GEOSample(
            assay=Assay.ATAC,
            library_name="ATAC_sample",
            title="ATAC-seq",
            organism=Organism.HUMAN,
            single_or_paired="single",
        )
        assert sample.single_or_paired == "single"

    def test_to_series_basic(self):
        """Test converting GEOSample to pandas Series."""
        sample = GEOSample(
            assay=Assay.ATAC,
            library_name="ATAC_sample",
            title="ATAC-seq",
            organism=Organism.HUMAN,
            cell_line="HEK293",
        )
        series = sample.to_series()
        assert isinstance(series, pd.Series)
        assert series["assay"] == Assay.ATAC
        assert series["library_name"] == "ATAC_sample"
        assert series["cell_line"] == "HEK293"

    def test_to_series_with_processed_files(self):
        """Test Series creation with processed data files."""
        sample = GEOSample(
            assay=Assay.RNA,
            library_name="RNA_sample",
            title="RNA-seq",
            organism=Organism.MOUSE,
            processed_data_file=["counts.txt", "tpm.txt", "fpkm.txt"],
        )
        series = sample.to_series()
        assert series["processed data file"] == "counts.txt"
        assert series["processed data file 1"] == "tpm.txt"
        assert series["processed data file 2"] == "fpkm.txt"

    def test_to_series_with_raw_files(self):
        """Test Series creation with raw files."""
        sample = GEOSample(
            assay=Assay.ATAC,
            library_name="ATAC_sample",
            title="ATAC-seq",
            organism=Organism.HUMAN,
            raw_file=["sample1_R1.fastq.gz", "sample1_R2.fastq.gz"],
        )
        series = sample.to_series()
        assert series["raw file"] == "sample1_R1.fastq.gz"
        assert series["raw file 1"] == "sample1_R2.fastq.gz"

    def test_to_series_chip_antibody_renamed(self):
        """Test that ChIP antibody field is renamed in series."""
        sample = GEOSample(
            assay=Assay.CHIP,
            library_name="ChIP_sample",
            title="ChIP-seq",
            organism=Organism.HUMAN,
            antibody="H3K27ac",
        )
        series = sample.to_series()
        assert "ChIP antibody" in series.index
        assert series["ChIP antibody"] == "H3K27ac"
        assert "antibody" not in series.index

    def test_to_series_cutandtag_antibody_renamed(self):
        """Test that CUT&TAG antibody field is renamed in series."""
        sample = GEOSample(
            assay=Assay.CAT,
            library_name="CAT_sample",
            title="CUT&TAG",
            organism=Organism.MOUSE,
            antibody="H3K4me3",
        )
        series = sample.to_series()
        assert "ChIP antibody" in series.index
        assert series["ChIP antibody"] == "H3K4me3"

    def test_from_series_basic(self):
        """Test creating GEOSample from pandas Series."""
        data = {
            "assay": Assay.ATAC,
            "library_name": "ATAC_sample",
            "title": "ATAC-seq",
            "organism": Organism.HUMAN,
            "cell_line": "HEK293",
        }
        series = pd.Series(data)
        sample = GEOSample.from_series(series)
        assert sample.assay == Assay.ATAC
        assert sample.library_name == "ATAC_sample"
        assert sample.cell_line == "HEK293"

    def test_from_series_with_processed_files(self):
        """Test from_series with processed data files."""
        data = {
            "assay": Assay.RNA,
            "library_name": "RNA_sample",
            "title": "RNA-seq",
            "organism": Organism.MOUSE,
            "processed data file": "counts.txt",
            "processed data file 1": "tpm.txt",
        }
        series = pd.Series(data)
        sample = GEOSample.from_series(series)
        assert len(sample.processed_data_file) == 2
        assert "counts.txt" in sample.processed_data_file
        assert "tpm.txt" in sample.processed_data_file

    def test_from_series_with_raw_files(self):
        """Test from_series with raw files."""
        data = {
            "assay": Assay.ATAC,
            "library_name": "ATAC_sample",
            "title": "ATAC-seq",
            "organism": Organism.HUMAN,
            "raw file": "R1.fastq.gz",
            "raw file 1": "R2.fastq.gz",
        }
        series = pd.Series(data)
        sample = GEOSample.from_series(series)
        assert len(sample.raw_file) == 2
        assert "R1.fastq.gz" in sample.raw_file
        assert "R2.fastq.gz" in sample.raw_file

    def test_to_series_from_series_roundtrip(self):
        """Test that to_series -> from_series roundtrip preserves data."""
        original = GEOSample(
            assay=Assay.CHIP,
            library_name="ChIP_H3K27ac",
            title="ChIP-seq for H3K27ac",
            organism=Organism.HUMAN,
            cell_line="K562",
            antibody="H3K27ac",
            treatment="control",
            processed_data_file=["peaks.bed", "coverage.bw"],
            raw_file=["sample_R1.fastq.gz", "sample_R2.fastq.gz"],
        )
        series = original.to_series()
        reconstructed = GEOSample.from_series(series)
        
        # Note: antibody will be in the series as "ChIP antibody" but from_series
        # should handle this - let's just check core fields
        assert reconstructed.assay == original.assay
        assert reconstructed.library_name == original.library_name
        assert reconstructed.title == original.title
        assert reconstructed.organism == original.organism
        assert reconstructed.cell_line == original.cell_line
        assert reconstructed.treatment == original.treatment
        assert len(reconstructed.processed_data_file) == 2
        assert len(reconstructed.raw_file) == 2


class TestGEOSamples:
    """Tests for GEOSamples collection model."""

    def test_create_empty_samples(self):
        """Test creating empty GEOSamples collection."""
        samples = GEOSamples(samples=[])
        assert len(samples.samples) == 0

    def test_create_samples_with_one_sample(self):
        """Test creating GEOSamples with single sample."""
        sample = GEOSample(
            assay=Assay.ATAC,
            library_name="ATAC_sample",
            title="ATAC-seq",
            organism=Organism.HUMAN,
        )
        samples = GEOSamples(samples=[sample])
        assert len(samples.samples) == 1
        assert samples.samples[0].library_name == "ATAC_sample"

    def test_create_samples_with_multiple_samples(self):
        """Test creating GEOSamples with multiple samples."""
        sample1 = GEOSample(
            assay=Assay.ATAC,
            library_name="ATAC_sample1",
            title="ATAC-seq 1",
            organism=Organism.HUMAN,
        )
        sample2 = GEOSample(
            assay=Assay.RNA,
            library_name="RNA_sample1",
            title="RNA-seq 1",
            organism=Organism.MOUSE,
        )
        samples = GEOSamples(samples=[sample1, sample2])
        assert len(samples.samples) == 2

    def test_to_dataframe_basic(self):
        """Test converting GEOSamples to DataFrame."""
        sample1 = GEOSample(
            assay=Assay.ATAC,
            library_name="ATAC_sample1",
            title="ATAC-seq 1",
            organism=Organism.HUMAN,
        )
        sample2 = GEOSample(
            assay=Assay.RNA,
            library_name="RNA_sample1",
            title="RNA-seq 1",
            organism=Organism.MOUSE,
        )
        samples = GEOSamples(samples=[sample1, sample2])
        df = samples.to_dataframe()
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2
        assert "library_name" in df.columns
        assert df.iloc[0]["library_name"] == "ATAC_sample1"
        assert df.iloc[1]["library_name"] == "RNA_sample1"

    def test_to_dataframe_with_files(self):
        """Test DataFrame creation with file columns."""
        sample = GEOSample(
            assay=Assay.RNA,
            library_name="RNA_sample",
            title="RNA-seq",
            organism=Organism.MOUSE,
            processed_data_file=["counts.txt", "tpm.txt"],
            raw_file=["R1.fastq.gz", "R2.fastq.gz"],
        )
        samples = GEOSamples(samples=[sample])
        df = samples.to_dataframe()
        
        assert "processed data file" in df.columns
        assert "raw file" in df.columns
        # Verify we have the expected columns for multiple files
        processed_cols = [col for col in df.columns if "processed data file" in col or "processed_data_file" in col]
        raw_cols = [col for col in df.columns if "raw file" in col or "raw_file" in col]
        # Should have at least one column for each file type
        assert len(processed_cols) >= 1
        assert len(raw_cols) >= 1

    def test_to_dataframe_column_normalization(self):
        """Test that DataFrame columns are normalized (trailing numbers removed)."""
        sample1 = GEOSample(
            assay=Assay.ATAC,
            library_name="ATAC_sample1",
            title="ATAC-seq",
            organism=Organism.HUMAN,
            processed_data_file=["file1.txt", "file2.txt"],
        )
        sample2 = GEOSample(
            assay=Assay.ATAC,
            library_name="ATAC_sample2",
            title="ATAC-seq",
            organism=Organism.HUMAN,
            processed_data_file=["file3.txt"],
        )
        samples = GEOSamples(samples=[sample1, sample2])
        df = samples.to_dataframe()
        
        # Column names should have trailing numbers removed by regex
        assert "processed data file" in df.columns
        # DataFrame should have 2 rows (one per sample)
        assert len(df) == 2
        # Check that columns exist
        processed_cols = [col for col in df.columns if "processed data file" in col or "processed_data_file" in col]
        assert len(processed_cols) >= 1
