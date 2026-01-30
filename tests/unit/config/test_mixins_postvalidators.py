"""Additional tests for post-validation behavior in config mixins."""

from pydantic import BaseModel

from seqnado.config.mixins import (
    CommonComputedFieldsMixin,
    PeakCallingMixin,
    SNPCallingMixin,
    MethylationMixin,
)


def test_create_bigwigs_false_for_methylation_mixin():
    """When model mixes in MethylationMixin, create_bigwigs must be False."""

    class TestModel(CommonComputedFieldsMixin, MethylationMixin, BaseModel):
        bigwigs: dict | None = None

    model = TestModel(bigwigs={"some": "config"})
    assert model.create_bigwigs is False


def test_create_bigwigs_respects_explicit_flag():
    """Explicitly provided create_bigwigs should be respected."""

    class TestModel(CommonComputedFieldsMixin, BaseModel):
        bigwigs: dict | None = None

    m1 = TestModel(bigwigs={"x": 1}, create_bigwigs=True)
    assert m1.create_bigwigs is True

    m2 = TestModel(bigwigs={"x": 1}, create_bigwigs=False)
    assert m2.create_bigwigs is False


def test_peak_calling_postvalidator_sets_flag():
    class TestModel(PeakCallingMixin, BaseModel):
        peak_calling: dict | None = None

    m_with = TestModel(peak_calling={"method": ["macs2"]})
    assert m_with.call_peaks is True

    m_without = TestModel(peak_calling=None)
    assert m_without.call_peaks is False


def test_snp_calling_postvalidator_sets_flag():
    class TestModel(SNPCallingMixin, BaseModel):
        snp_calling: dict | None = None

    m_with = TestModel(snp_calling={"method": "bcftools"})
    assert m_with.call_snps is True

    m_without = TestModel(snp_calling=None)
    assert m_without.call_snps is False


def test_methylation_postvalidator_sets_flag():
    class TestModel(MethylationMixin, BaseModel):
        methylation: dict | None = None

    m_with = TestModel(methylation={"method": "taps"})
    assert m_with.call_methylation is True

    m_without = TestModel(methylation=None)
    assert m_without.call_methylation is False
