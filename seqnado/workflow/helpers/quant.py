"""Helper functions for quantification workflows."""


def get_bams_to_count(CONFIG, SAMPLE_NAMES, OUTPUT_DIR):
    """
    Get list of BAM files to use for counting.

    Returns spike-in filtered BAMs if spike-in is configured, otherwise returns
    regular aligned BAMs.

    Args:
        CONFIG: The configuration object.
        SAMPLE_NAMES: List of sample names.
        OUTPUT_DIR: The output directory path.

    Returns:
        list: List of BAM file paths to count.
    """
    if CONFIG.assay_config.has_spikein:
        return [
            OUTPUT_DIR + f"/aligned/spikein/filtered/{sample}.bam"
            for sample in SAMPLE_NAMES
        ]
    return [OUTPUT_DIR + f"/aligned/{sample}.bam" for sample in SAMPLE_NAMES]
