from itertools import chain

def get_mcc_bigwig_files(wildcards):
    """
    Get MCC bigwig files for a given sample.

    Args:
        wildcards: Snakemake wildcards object containing 'sample'.
    Returns:
        list: List of bigwig file paths for the sample.
    """
    return list(chain.from_iterable(rules.confirm_bigwigs_generated.input))

def get_mcc_peak_files(wildcards):
    """
    Get MCC peak files for a given sample.

    Args:
        wildcards: Snakemake wildcards object containing 'sample'.
    Returns:
        list: List of peak file paths for the sample.
    """
    return list(chain.from_iterable(rules.confirm_peaks_generated.input))

def get_hub_input_files(wildcards, OUTPUT, ASSAY):
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
        input_files.extend(get_mcc_bigwig_files(wildcards))
        input_files.extend(get_mcc_peak_files(wildcards))
    
    return input_files
