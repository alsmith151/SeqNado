"""Tests for seqnado.config.mixins module."""

import pytest
from pydantic import BaseModel

from seqnado.config.mixins import (
    CommonComputedFieldsMixin,
    PeakCallingMixin,
    SNPCallingMixin,
    MethylationMixin,
    PathValidatorMixin,
)


class TestCommonComputedFieldsMixin:
    """Tests for CommonComputedFieldsMixin computed fields."""

    def test_create_bigwigs_when_config_present(self):
        """Test create_bigwigs returns True when bigwigs config is set."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            bigwigs: dict | None = None
            create_geo_submission_files: bool = False
        
        model = TestModel(bigwigs={"some": "config"})
        assert model.create_bigwigs is True

    def test_create_bigwigs_when_config_absent(self):
        """Test create_bigwigs returns False when bigwigs config is None."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            bigwigs: dict | None = None
            create_geo_submission_files: bool = False
        
        model = TestModel(bigwigs=None)
        assert model.create_bigwigs is False

    def test_plot_with_plotnado_when_config_present(self):
        """Test plot_with_plotnado returns True when plotting config is set."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            plotting: dict | None = None
            create_geo_submission_files: bool = False
        
        model = TestModel(plotting={"some": "config"})
        assert model.plot_with_plotnado is True

    def test_plot_with_plotnado_when_config_absent(self):
        """Test plot_with_plotnado returns False when plotting config is None."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            plotting: dict | None = None
            create_geo_submission_files: bool = False
        
        model = TestModel(plotting=None)
        assert model.plot_with_plotnado is False

    def test_create_dataset_when_config_present(self):
        """Test create_dataset returns True when dataset_for_ml is set."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            dataset_for_ml: dict | None = None
            create_geo_submission_files: bool = False
        
        model = TestModel(dataset_for_ml={"binsize": 1000})
        assert model.create_dataset is True

    def test_create_dataset_when_config_absent(self):
        """Test create_dataset returns False when dataset_for_ml is None."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            dataset_for_ml: dict | None = None
            create_geo_submission_files: bool = False
        
        model = TestModel(dataset_for_ml=None)
        assert model.create_dataset is False

    def test_create_ucsc_hub_when_config_present(self):
        """Test create_ucsc_hub returns True when ucsc_hub is set."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            ucsc_hub: dict | None = None
            create_geo_submission_files: bool = False
        
        model = TestModel(ucsc_hub={"directory": "hub/"})
        assert model.create_ucsc_hub is True

    def test_create_ucsc_hub_when_config_absent(self):
        """Test create_ucsc_hub returns False when ucsc_hub is None."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            ucsc_hub: dict | None = None
            create_geo_submission_files: bool = False
        
        model = TestModel(ucsc_hub=None)
        assert model.create_ucsc_hub is False

    def test_has_spikein_when_config_present(self):
        """Test has_spikein returns True when spikein is set."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            spikein: dict | None = None
            create_geo_submission_files: bool = False
        
        model = TestModel(spikein={"genome": "dm6"})
        assert model.has_spikein is True

    def test_has_spikein_when_config_absent(self):
        """Test has_spikein returns False when spikein is None."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            spikein: dict | None = None
            create_geo_submission_files: bool = False
        
        model = TestModel(spikein=None)
        assert model.has_spikein is False

    def test_validate_geo_submission_files_true(self):
        """Test geo submission files validator with True value."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            create_geo_submission_files: bool = False
        
        model = TestModel(create_geo_submission_files=True)
        assert model.create_geo_submission_files is True

    def test_validate_geo_submission_files_false(self):
        """Test geo submission files validator with False value."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            create_geo_submission_files: bool = False
        
        model = TestModel(create_geo_submission_files=False)
        assert model.create_geo_submission_files is False


class TestPeakCallingMixin:
    """Tests for PeakCallingMixin."""

    def test_call_peaks_when_config_present(self):
        """Test call_peaks returns True when peak_calling config is set."""
        
        class TestModel(PeakCallingMixin, BaseModel):
            peak_calling: dict | None = None
        
        model = TestModel(peak_calling={"method": ["macs2"]})
        assert model.call_peaks is True

    def test_call_peaks_when_config_absent(self):
        """Test call_peaks returns False when peak_calling config is None."""
        
        class TestModel(PeakCallingMixin, BaseModel):
            peak_calling: dict | None = None
        
        model = TestModel(peak_calling=None)
        assert model.call_peaks is False


class TestSNPCallingMixin:
    """Tests for SNPCallingMixin."""

    def test_call_snps_when_config_present(self):
        """Test call_snps returns True when snp_calling config is set."""
        
        class TestModel(SNPCallingMixin, BaseModel):
            snp_calling: dict | None = None
        
        model = TestModel(snp_calling={"method": "bcftools"})
        assert model.call_snps is True

    def test_call_snps_when_config_absent(self):
        """Test call_snps returns False when snp_calling config is None."""
        
        class TestModel(SNPCallingMixin, BaseModel):
            snp_calling: dict | None = None
        
        model = TestModel(snp_calling=None)
        assert model.call_snps is False


class TestMethylationMixin:
    """Tests for MethylationMixin."""

    def test_call_methylation_when_config_present(self):
        """Test call_methylation returns True when methylation config is set."""
        
        class TestModel(MethylationMixin, BaseModel):
            methylation: dict | None = None
        
        model = TestModel(methylation={"method": "taps"})
        assert model.call_methylation is True

    def test_call_methylation_when_config_absent(self):
        """Test call_methylation returns False when methylation config is None."""
        
        class TestModel(MethylationMixin, BaseModel):
            methylation: dict | None = None
        
        model = TestModel(methylation=None)
        assert model.call_methylation is False


class TestPathValidatorMixin:
    """Tests for PathValidatorMixin."""

    def test_validate_path_exists_with_existing_path(self, tmp_path):
        """Test path validation with an existing path."""
        test_file = tmp_path / "test.txt"
        test_file.touch()
        
        result = PathValidatorMixin.validate_path_exists(test_file, "test_path")
        assert result == test_file

    def test_validate_path_exists_with_existing_path_string(self, tmp_path):
        """Test path validation with an existing path as string."""
        test_file = tmp_path / "test.txt"
        test_file.touch()
        
        result = PathValidatorMixin.validate_path_exists(str(test_file), "test_path")
        assert result == str(test_file)

    def test_validate_path_exists_with_nonexistent_path(self):
        """Test path validation with a non-existent path raises error."""
        with pytest.raises(ValueError, match="test_path /nonexistent/path does not exist"):
            PathValidatorMixin.validate_path_exists("/nonexistent/path", "test_path")

    def test_validate_path_exists_with_none(self):
        """Test path validation with None returns None."""
        result = PathValidatorMixin.validate_path_exists(None, "test_path")
        assert result is None

    def test_validate_required_when_condition_true_with_value(self):
        """Test required validation when condition is met and value provided."""
        result = PathValidatorMixin.validate_required_when(
            "some_value", condition=True, field_name="test_field"
        )
        assert result == "some_value"

    def test_validate_required_when_condition_true_without_value(self):
        """Test required validation when condition is met but value is missing."""
        with pytest.raises(ValueError, match="test_field must be provided"):
            PathValidatorMixin.validate_required_when(
                None, condition=True, field_name="test_field"
            )

    def test_validate_required_when_condition_false(self):
        """Test required validation when condition is not met."""
        result = PathValidatorMixin.validate_required_when(
            None, condition=False, field_name="test_field"
        )
        assert result is None


class TestComputedFieldOverrides:
    """Tests for computed field override functionality."""

    def test_create_bigwigs_explicit_override_true(self):
        """Test create_bigwigs can be explicitly set to True (with config present)."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            bigwigs: dict | None = None
            create_geo_submission_files: bool = False
        
        # Explicitly set create_bigwigs=True with bigwigs config
        model = TestModel(bigwigs={"some": "config"}, create_bigwigs=True)
        assert model.create_bigwigs is True

    def test_create_bigwigs_explicit_override_false(self):
        """Test create_bigwigs can be explicitly set to False even when config is present."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            bigwigs: dict | None = None
            create_geo_submission_files: bool = False
        
        # Explicitly set create_bigwigs=False despite bigwigs config
        model = TestModel(bigwigs={"some": "config"}, create_bigwigs=False)
        assert model.create_bigwigs is False

    def test_create_bigwigs_override_validation_error(self):
        """Test create_bigwigs=True without bigwigs config raises validation error."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            bigwigs: dict | None = None
            create_geo_submission_files: bool = False
        
        # Should raise error when create_bigwigs=True but no bigwigs config
        with pytest.raises(ValueError, match="create_bigwigs=True requires 'bigwigs' configuration"):
            TestModel(bigwigs=None, create_bigwigs=True)

    def test_plot_with_plotnado_explicit_override(self):
        """Test plot_with_plotnado can be explicitly overridden."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            plotting: dict | None = None
            create_geo_submission_files: bool = False
        
        # Override to False despite config present
        model = TestModel(plotting={"some": "config"}, plot_with_plotnado=False)
        assert model.plot_with_plotnado is False

    def test_plot_with_plotnado_override_validation_error(self):
        """Test plot_with_plotnado=True without plotting config raises error."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            plotting: dict | None = None
            create_geo_submission_files: bool = False
        
        with pytest.raises(ValueError, match="plot_with_plotnado=True requires 'plotting' configuration"):
            TestModel(plotting=None, plot_with_plotnado=True)

    def test_create_dataset_override_validation_error(self):
        """Test create_dataset=True without dataset_for_ml config raises error."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            dataset_for_ml: dict | None = None
            create_geo_submission_files: bool = False
        
        with pytest.raises(ValueError, match="create_dataset=True requires 'dataset_for_ml' configuration"):
            TestModel(dataset_for_ml=None, create_dataset=True)

    def test_create_ucsc_hub_override_validation_error(self):
        """Test create_ucsc_hub=True without ucsc_hub config raises error."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            ucsc_hub: dict | None = None
            create_geo_submission_files: bool = False
        
        with pytest.raises(ValueError, match="create_ucsc_hub=True requires 'ucsc_hub' configuration"):
            TestModel(ucsc_hub=None, create_ucsc_hub=True)

    def test_has_spikein_override_validation_error(self):
        """Test has_spikein=True without spikein config raises error."""
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            spikein: dict | None = None
            create_geo_submission_files: bool = False
        
        with pytest.raises(ValueError, match="has_spikein=True requires 'spikein' configuration"):
            TestModel(spikein=None, has_spikein=True)

    def test_call_peaks_override_validation_error(self):
        """Test call_peaks=True without peak_calling config raises error."""
        
        class TestModel(PeakCallingMixin, BaseModel):
            peak_calling: dict | None = None
        
        with pytest.raises(ValueError, match="call_peaks=True requires 'peak_calling' configuration"):
            TestModel(peak_calling=None, call_peaks=True)

    def test_call_snps_override_validation_error(self):
        """Test call_snps=True without snp_calling config raises error."""
        
        class TestModel(SNPCallingMixin, BaseModel):
            snp_calling: dict | None = None
        
        with pytest.raises(ValueError, match="call_snps=True requires 'snp_calling' configuration"):
            TestModel(snp_calling=None, call_snps=True)

    def test_call_methylation_override_validation_error(self):
        """Test call_methylation=True without methylation config raises error."""
        
        class TestModel(MethylationMixin, BaseModel):
            methylation: dict | None = None
        
        with pytest.raises(ValueError, match="call_methylation=True requires 'methylation' configuration"):
            TestModel(methylation=None, call_methylation=True)

    def test_yaml_roundtrip_with_overrides(self, tmp_path):
        """Test that YAML serialization and deserialization preserves explicit overrides."""
        from seqnado.config.yaml_serializer import dump_config_to_yaml, load_config_from_yaml
        
        class TestModel(CommonComputedFieldsMixin, BaseModel):
            bigwigs: dict | None = None
            plotting: dict | None = None
            create_geo_submission_files: bool = False
        
        # Create model with explicit override
        original = TestModel(
            bigwigs={"some": "config"}, 
            plotting={"other": "config"},
            create_bigwigs=False,  # Override to False despite config
            plot_with_plotnado=True  # Explicitly True (matches computed)
        )
        
        # Serialize to YAML
        yaml_file = tmp_path / "test_config.yaml"
        dump_config_to_yaml(original, stream=yaml_file)
        
        # Load back from YAML
        data = load_config_from_yaml(yaml_file)
        restored = TestModel(**data)
        
        # Verify overrides are preserved
        assert restored.create_bigwigs is False  # Overridden value
        assert restored.plot_with_plotnado is True  # Explicit value
