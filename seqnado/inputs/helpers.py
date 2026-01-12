from pathlib import Path
from seqnado import Assay
import pandas as pd
from .fastq import FastqCollection, FastqCollectionForIP
from .bam import BamCollection
from .bigwigs import BigWigCollection
from .interfaces import CollectionLike


def get_sample_collection(assay: Assay, path: Path, config=None) -> CollectionLike:
    """
    Select the appropriate sample collection class based on the contents of the metadata file.

    Args:
        assay (Assay): The assay type (e.g., 'ChIP-seq', 'RNA-seq').
        path (Path): Path to the metadata file
        config: Optional config object with has_spikein property to determine if deseq2 validation is required

    Returns:
        Type[Union[FastqCollection, FastqCollectionForIP, BamCollection, BigWigCollection]]: The appropriate collection class.
    Raises:
        RuntimeError: If the type of collection cannot be determined from the metadata.
    """
    # Step 1: Read the metadata file
    if not path.exists():
        raise FileNotFoundError(f"Metadata file not found: {path}")
    
    df = pd.read_csv(path)
    columns = df.columns.str.lower()

    # Step 2: Determine the type of collection based on columns
    if any(col in columns for col in ['r1', 'r2']):
        if 'ip' in columns:
            type_of_collection = FastqCollectionForIP
        else:
            type_of_collection = FastqCollection
    elif 'bam' in columns:
        type_of_collection = BamCollection

    # For bigWig files, we check for 'bigwig' or strand-specific columns
    elif any(col in columns for col in ['bigwig', 'bigwig_plus', 'bigwig_minus']):
        type_of_collection = BigWigCollection
    
    else:
        raise RuntimeError("Could not determine the type of collection from metadata columns.")

    # Determine if deseq2 validation should be enforced (only for RNA with spike-in)
    validate_deseq2 = getattr(config, 'has_spikein', False) if config else False
    
    return type_of_collection.from_csv(assay=assay, path=path, validate_deseq2=validate_deseq2, assay_for_validation=assay)
