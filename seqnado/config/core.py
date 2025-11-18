from pathlib import Path
from typing import Union, Optional
from enum import Enum
from pydantic import BaseModel, computed_field, field_validator, Field, model_validator
from seqnado import Assay
from .configs import (
    BigwigConfig,
    PlottingConfig,
    UCSCHubConfig,
    MLDatasetConfig,
    RNAQuantificationConfig,
    PCRDuplicatesConfig,
    QCConfig,
    GenomeConfig,
    ProjectConfig,
    PeakCallingConfig,
    SNPCallingConfig,
    SpikeInConfig,
    MethylationConfig,
    MCCConfig,
)
from .mixins import (
    CommonComputedFieldsMixin,
    PeakCallingMixin,
    SNPCallingMixin,
    MethylationMixin,
)
from .third_party_tools import ThirdPartyToolsConfig


class BaseAssayConfig(BaseModel, CommonComputedFieldsMixin):
    """Base configuration for all assays."""

    bigwigs: BigwigConfig | None = None
    plotting: PlottingConfig | None = None
    ucsc_hub: UCSCHubConfig | None = None
    dataset_for_ml: MLDatasetConfig | None = None

    # Boolean flags for optional features
    create_geo_submission_files: bool = False


class ATACAssayConfig(BaseAssayConfig, PeakCallingMixin):
    """Configuration specific to ATAC-seq assays."""

    tn5_shift: bool = False
    peak_calling: PeakCallingConfig | None = None
    create_heatmaps: bool = False


class ChIPAssayConfig(BaseAssayConfig, PeakCallingMixin):
    """Configuration specific to ChIP-seq assays."""

    spikein: SpikeInConfig | None = None
    peak_calling: PeakCallingConfig | None = None
    create_heatmaps: bool = False


class CATAssayConfig(BaseAssayConfig, PeakCallingMixin):
    """Configuration specific to CAT-seq assays."""

    spikein: SpikeInConfig | None = None
    tn5_shift: bool = False
    peak_calling: PeakCallingConfig | None = None
    create_heatmaps: bool = False


class RNAAssayConfig(BaseAssayConfig):
    """Configuration specific to RNA-seq assays."""

    rna_quantification: RNAQuantificationConfig | None = None
    create_heatmaps: bool = False


class SNPAssayConfig(BaseAssayConfig, SNPCallingMixin):
    """Configuration specific to SNP calling assays."""

    snp_calling: SNPCallingConfig | None = None
    ucsc_hub: None


class MCCAssayConfig(BaseAssayConfig):
    """Configuration specific to MCC (Capture-C) assays."""

    mcc: MCCConfig | None = None


class MethylationAssayConfig(BaseAssayConfig, MethylationMixin):
    """Configuration specific to methylation assays."""

    methylation: MethylationConfig | None = None
    ucsc_hub: None


class CRISPRAssayConfig(BaseAssayConfig):
    """Configuration specific to CRISPR assays."""

    # CRISPR-specific options can be added here
    ucsc_hub: None


# Union type for all assay-specific configurations
AssaySpecificConfig = Union[
    ATACAssayConfig,
    ChIPAssayConfig,
    CATAssayConfig,
    RNAAssayConfig,
    SNPAssayConfig,
    MCCAssayConfig,
    MethylationAssayConfig,
    CRISPRAssayConfig,
]

# Class constant for assay config mapping to reduce redundancy
ASSAY_CONFIG_MAP = {
    Assay.ATAC: ATACAssayConfig,
    Assay.CHIP: ChIPAssayConfig,
    Assay.CAT: CATAssayConfig,
    Assay.RNA: RNAAssayConfig,
    Assay.SNP: SNPAssayConfig,
    Assay.MCC: MCCAssayConfig,
    Assay.METH: MethylationAssayConfig,
    Assay.CRISPR: CRISPRAssayConfig,
}


class AssayConfig(Enum):
    ATAC = ATACAssayConfig
    CHIP = ChIPAssayConfig
    CAT = CATAssayConfig
    RNA = RNAAssayConfig
    SNP = SNPAssayConfig
    MCC = MCCAssayConfig
    METH = MethylationAssayConfig
    CRISPR = CRISPRAssayConfig


class SeqnadoConfig(BaseModel):
    """Configuration for the SeqNado workflow."""

    assay: Assay
    project: ProjectConfig
    genome: GenomeConfig
    metadata: Path
    qc: QCConfig = QCConfig()
    pcr_duplicates: PCRDuplicatesConfig = PCRDuplicatesConfig()
    remove_blacklist: bool = False
    assay_config: AssaySpecificConfig | None = None
    third_party_tools: ThirdPartyToolsConfig | None = Field(None, description="Configuration for third-party tools.")

    # If no third_party_tools config is provided, use defaults
    @model_validator(mode="before")
    def set_default_third_party_tools(cls, values):
        if "third_party_tools" not in values or values["third_party_tools"] is None:
            values["third_party_tools"] = ThirdPartyToolsConfig.for_assay(values.get("assay"))
        return values

    @classmethod
    def from_yaml(cls, path: Path) -> "SeqnadoConfig":
        """Load configuration from a YAML file."""
        import yaml

        with open(path, "r") as f:
            data = yaml.safe_load(f)

        return cls(**data)

    @computed_field
    @property
    def organism(self) -> str:
        """Return the organism (string) from the genome configuration."""
        return self.genome.organism

    @computed_field
    @property
    def shift_for_tn5_insertion(self) -> bool:
        """Return the Tn5 shift configuration for the specified assay."""
        return hasattr(self.assay_config, "tn5_shift") and self.assay_config.tn5_shift
    
    @computed_field
    @property
    def mcc_viewpoints(self) -> str:
        """Return the MCC viewpoints file path."""
        if self.assay_config and hasattr(self.assay_config, "mcc"):
            return str(self.assay_config.mcc.viewpoints)
        return ""

    @field_validator("remove_blacklist")
    def validate_remove_blacklist(cls, v):
        """Can only be set to True if genome blacklist is provided."""
        if v and not cls.genome.blacklist:
            raise ValueError(
                "remove_blacklist can only be True if genome blacklist is provided."
            )
        return v

    @field_validator("assay_config", mode="before")
    def validate_assay_config_matches_assay(cls, v, info):
        """Ensure the assay_config type matches the specified assay."""
        if v is None:
            return v

        assay = info.data.get("assay")
        if assay is None:
            return v

        expected_config_class = ASSAY_CONFIG_MAP.get(assay)
        if expected_config_class and not isinstance(v, expected_config_class):
            if isinstance(v, dict):
                # Try to create the appropriate config from dict
                return expected_config_class(**v)
            else:
                raise ValueError(
                    f"assay_config must be of type {expected_config_class.__name__} for assay {assay.value}"
                )

        return v

    @classmethod
    def create_assay_config(cls, assay: Assay, **kwargs) -> AssaySpecificConfig:
        """Create the appropriate assay config for the given assay type."""
        config_class = cls._ASSAY_CONFIG_MAP.get(assay)
        if config_class is None:
            raise ValueError(
                f"No configuration class available for assay {assay.value}"
            )

        return config_class(**kwargs)
    
    @classmethod
    def from_yaml(cls, path: Path) -> "SeqnadoConfig":
        """Load configuration from a YAML file."""
        import yaml

        with open(path, "r") as f:
            data = yaml.safe_load(f)

        return cls(**data)
