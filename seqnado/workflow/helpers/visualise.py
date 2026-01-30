from itertools import chain
from seqnado import Assay

def get_mcc_bigwig_files(wildcards, rules):
    """
    Get MCC bigwig files for a given sample.

    Args:
        wildcards: Snakemake wildcards object containing 'sample'.
    Returns:
        list: List of bigwig file paths for the sample.
    """
    return rules.confirm_bigwigs_generated.input

def get_mcc_peak_files(wildcards, rules):
    """
    Get MCC peak files for a given sample.

    Args:
        wildcards: Snakemake wildcards object containing 'sample'.
    Returns:
        list: List of peak file paths for the sample.
    """
    return rules.confirm_peaks_generated.input

def get_hub_input_files(wildcards, OUTPUT, ASSAY, rules):
    """
    Get all input files for UCSC hub generation, including MCC bigwigs and peaks if applicable.

    Args:
        wildcards: Snakemake wildcards object.

    Returns:
        list: List of file paths to include in the UCSC hub.
    """
    input_files = []
    input_files.extend(OUTPUT.select_files(suffix=".bigwig", exclude=["/geo_submission/"]))
    input_files.extend(OUTPUT.bigbed_files)
    
    if ASSAY == Assay.MCC:
        input_files.extend(get_mcc_bigwig_files(wildcards, rules))
        input_files.extend(get_mcc_peak_files(wildcards, rules))
    
    return input_files
