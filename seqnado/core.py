from enum import Enum


# =============================================================================
# CONSTANTS
# =============================================================================
NONE_VALUES = [None, "None", "none", "null", "Null", "NULL", ".", "", "NA"]

# ============================================================
# ENUMS
# ============================================================


class Assay(Enum):
    """Supported sequencing assay types."""

    RNA = "rna"
    ATAC = "atac"
    SNP = "snp"
    CHIP = "chip"
    CAT = "cat"
    METH = "meth"
    MCC = "mcc"
    CRISPR = "crispr"

    @classmethod
    def non_ip_assays(cls):
        """Return assays that don't require IP (immunoprecipitation)."""
        ip_assays = {cls.CHIP, cls.CAT}
        return [assay for assay in cls if assay not in ip_assays]

    @classmethod
    def ip_assays(cls):
        """Return assays that require IP (immunoprecipitation)."""
        return [cls.CHIP, cls.CAT]


class PileupMethod(Enum):
    """Methods for creating pileup files."""

    DEEPTOOLS = "deeptools"
    HOMER = "homer"
    BAMNADO = "bamnado"


class ScaleMethod(Enum):
    """Methods for scaling genomic data."""

    UNSCALED = "unscaled"
    CSAW = "csaw"
    CPM = "cpm"
    RPKM = "rpkm"
    SPIKEIN = "spikein"
    MERGED = "merged"


class PeakCallingMethod(Enum):
    """Methods for calling peaks."""

    MACS2 = "macs2"
    MACS3 = "macs3"
    HOMER = "homer"
    LANCEOTRON = "lanceotron"
    SEACR = "seacr"
    LANCEOTRON_MCC = "lanceotron-mcc"


class PCRDuplicateHandling(Enum):
    """Methods for handling PCR duplicates."""

    REMOVE = "remove"
    MARK = "mark"
    NONE = "dont_remove"


class PCRDuplicateTool(Enum):
    PICARD = "picard"
    SAMTOOLS = "samtools"



class SpikeInMethod(Enum):
    """Methods for spike-in normalization."""

    ORLANDO = "orlando"
    WITH_INPUT = "with_input"


class SNPCallingMethod(Enum):
    """Methods for SNP calling."""

    BCFTOOLS = "bcftools"
    DEEPVARIANT = "deepvariant"


class RNAQuantificationMethod(Enum):
    """Methods for RNA quantification."""

    FEATURE_COUNTS = "feature_counts"
    SALMON = "salmon"


class MethylationMethod(Enum):
    """Methods for methylation calling."""

    TAPS = "taps"
    BISULFITE = "bisulfite"