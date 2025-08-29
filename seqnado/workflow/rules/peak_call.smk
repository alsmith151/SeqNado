from seqnado import FileType, DataScalingTechnique
from seqnado.inputs import FastqCollection, FastqCollectionForIP
from seqnado.config.third_party_tools import CommandLineArguments 


def get_control_file(wildcards, file_type: FileType):
    """
    Get the control file for a given treatment sample.

    Args:
        wildcards: Snakemake wildcards object containing 'sample' and 'treatment'.
        file_type (FileType): The type of file to retrieve (e.g., FileType.BAM, FileType.BIGWIG, FileType.TAG).

    Returns:
        str or None: The path to the control file if it exists, otherwise None.
    """

    if isinstance(INPUT_FILES, FastqCollection):
        return "NO_CONTROL_PRESENT_FOR_FILE" # Dummy value to prevent the rule from ever being able to run
    elif isinstance(INPUT_FILES, FastqCollectionForIP):
        search_term = f"{wildcards.sample}_{wildcards.treatment}"
        control_name = INPUT_FILES.get_control_performed(search_term)
        if not control_name:
            return "NO_CONTROL_PRESENT_FOR_FILE"

    
    # Match control with correct filetype
    match file_type:
        case FileType.BAM:
            return f"seqnado_output/aligned/{control_name}.bam"
        case FileType.BIGWIG:
            return f"seqnado_output/bigwigs/deeptools/{DataScalingTechnique.UNSCALED.value}/{control_name}.bigwig"
        case FileType.TAG:
            return f"seqnado_output/tag_dirs/{control_name}.tag"
        case _:
            raise ValueError(f"Unsupported file type '{file_type}' for control file retrieval")


def correct_macs_options(wildcards, options: CommandLineArguments):

    # Correct the options based on whether the input is paired or not
    sample_id = wildcards.sample
    treatment = wildcards.treatment if hasattr(wildcards, "treatment") else None
    search_term = f"{sample_id}_{treatment}" if treatment else sample_id
    is_paired = INPUT_FILES.is_paired(search_term)
    if is_paired:
        options = CommandLineArguments(value=options.value, exclude={"-f"})
    return options


rule macs2_with_input:
    input:
        treatment="seqnado_output/aligned/{sample}_{treatment}.bam",
        control=lambda wc: get_control_file(wc, file_type=FileType.BAM),
    output:
        peaks="seqnado_output/peaks/macs/{sample}_{treatment}.bed",
    params:
        options=lambda wc: format_macs_options(wc, config["macs"]["callpeak"]),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/macs/{sample}_{treatment}.log",
    container:
        "docker://quay.io/biocontainers/macs2:2.1.1.20160309--py27r3.3.1_1"
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} -n {params.basename} {params.options} > {log} 2>&1 &&
        cat {params.raw} | grep -v '^#' | grep -vE '^chr\\s+start\\s+end.*' | grep -v '^$' | cut -f 1-3 > {output.peaks}
        """


rule macs2_no_input:
    input:
        treatment="seqnado_output/aligned/{sample}_{treatment}.bam",
        control=None,
    output:
        peaks="seqnado_output/peaks/macs/{sample}_{treatment}.bed",
    params:
        options=lambda wc: str(format_macs_options(wc, CONFIG.third_party_tools.macs.callpeak.command_line_arguments)),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/macs/{sample}_{treatment}.log",
    container:
        "docker://quay.io/biocontainers/macs2:2.1.1.20160309--py27r3.3.1_1"
    shell:
        """
        macs2 callpeak -t {input.treatment} -n {params.basename} {params.options} > {log} 2>&1 &&
        cat {params.raw} | grep -v '^#' | grep -vE '^chr\\s+start\\s+end.*' | grep -v '^$' | cut -f 1-3 > {output.peaks}
        """


ruleorder: macs2_with_input > macs2_no_input

