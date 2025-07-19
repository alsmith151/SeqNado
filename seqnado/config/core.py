
from enum import Enum

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
    NONE = "none"

class PCRDuplicateTool(Enum):
    PICARD = "picard"
    SAMTOOLS = "samtools"
    