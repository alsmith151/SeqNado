from seqnado import FileType, DataScalingTechnique, Assay
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
        return "NO_CONTROL_PRESENT_FOR_FILE" # Dummy value to prevent the rule from ever being able to run if not an IP
    elif isinstance(INPUT_FILES, FastqCollectionForIP):
        search_term = wildcards.sample_id
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
    sample_id = wildcards.sample_id
    # treatment = wildcards.treatment if hasattr(wildcards, "treatment") else None
    # search_term = f"{sample_id}_{treatment}" if treatment else sample_id
    is_paired = INPUT_FILES.is_paired_end(sample_id)
    if not is_paired:
        options = CommandLineArguments(value=options.value, exclude={"-f"})
    return options


def get_lanceotron_call_peaks_threshold(wildcards):
    options = config["lanceotron"]["callpeak"]
    threshold_pattern = re.compile(r"\-c\s+(\d+.?\d*)")
    threshold = threshold_pattern.search(options).group(1)
    return threshold




rule macs2_with_input:
    input:
        treatment="seqnado_output/aligned/{sample_id}.bam",
        control=lambda wc: get_control_file(wc, file_type=FileType.BAM),
    output:
        peaks="seqnado_output/peaks/macs2/{sample_id}.bed",
    params:
        options=lambda wc: str(correct_macs_options(wc, CONFIG.third_party_tools.macs.call_peak.command_line_arguments)),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/macs2/{sample_id}.log",
    container:
        "docker://quay.io/biocontainers/macs2:2.1.1.20160309--py27r3.3.1_1"
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} -n {params.basename} {params.options} > {log} 2>&1 &&
        cat {params.raw} | grep -v '^#' | grep -vE '^chr\\s+start\\s+end.*' | grep -v '^$' | cut -f 1-3 > {output.peaks}
        """


rule macs2_no_input:
    input:
        treatment="seqnado_output/aligned/{sample_id}.bam",
    output:
        peaks="seqnado_output/peaks/macs2/{sample_id}.bed",
    params:
        options=lambda wc: str(correct_macs_options(wc, CONFIG.third_party_tools.macs.call_peak.command_line_arguments)),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/macs2/{sample_id}.log",
    container:
        "docker://quay.io/biocontainers/macs2:2.1.1.20160309--py27r3.3.1_1"
    shell:
        """
        macs2 callpeak -t {input.treatment} -n {params.basename} {params.options} > {log} 2>&1 &&
        cat {params.raw} | grep -v '^#' | grep -vE '^chr\\s+start\\s+end.*' | grep -v '^$' | cut -f 1-3 > {output.peaks}
        """


ruleorder: macs2_with_input > macs2_no_input



rule homer_with_input:
    input:
        treatment="seqnado_output/tag_dirs/{sample_id}",
        control=lambda wc: get_control_file(wc, file_type=FileType.TAG),
    output:
        peaks="seqnado_output/peaks/homer/{sample_id}.bed",
    log:
        "seqnado_output/logs/homer/{sample_id}.log",
    params:
        options=str(CONFIG.third_party_tools.homer.find_peaks.command_line_arguments),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp  -i {input.control} > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule homer_no_input:
    input:
        treatment="seqnado_output/tag_dirs/{sample_id}",
    output:
        peaks="seqnado_output/peaks/homer/{sample_id}.bed",
    log:
        "seqnado_output/logs/homer/{sample_id}.log",
    params:
        options=str(CONFIG.third_party_tools.homer.find_peaks.command_line_arguments),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """

ruleorder: homer_with_input > homer_no_input


rule lanceotron_with_input:
    input:
        treatment="seqnado_output/bigwigs/deeptools/unscaled/{sample_id}.bigWig",
        control=lambda wc: get_control_file(wc, file_type=FileType.BIGWIG),
    output:
        peaks="seqnado_output/peaks/lanceotron/{sample_id}.bed",
        ltron_peaks=temp("seqnado_output/peaks/lanceotron/{sample_id}_L-tron.bed"),
    log:
        "seqnado_output/logs/lanceotron/{sample_id}.log",
    params:
        threshold=lambda wildcards: get_lanceotron_call_peaks_threshold(wildcards),
        options=str(CONFIG.third_party_tools.lanceotron.call_peaks.command_line_arguments),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    container:
        "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest"
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=12, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    shell:"""
    lanceotron callPeaksInput {input.treatment} -i {input.control} -f {params.outdir} --skipheader > {log} 2>&1 &&
    cat {output.ltron_peaks} | awk 'BEGIN{{OFS="\\t"}} $4 >= {params.threshold} {{print $1, $2, $3}}' > {output.peaks} 
    """


rule lanceotron_no_input:
    input:
        treatment="seqnado_output/bigwigs/deeptools/unscaled/{sample_id}.bigWig",
    output:
        peaks="seqnado_output/peaks/lanceotron/{sample_id}.bed",
        ltron_peaks=temp("seqnado_output/peaks/lanceotron/{sample_id}_L-tron.bed"),
    log:
        "seqnado_output/logs/lanceotron/{sample_id}.log",
    params:
        options=str(CONFIG.third_party_tools.lanceotron.call_peaks.command_line_arguments),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    container:
        "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=12, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    shell:"""
    lanceotron callPeaks {input.treatment} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
    cat {output.ltron_peaks} | cut -f 1-3 > {output.peaks}
    """

ruleorder: lanceotron_with_input > lanceotron_no_input


if ASSAY == Assay.CAT:

    rule seacr:
        input:
            treatment="seqnado_output/bedgraphs/{sample_id}.bedGraph",
        output:
            peaks="seqnado_output/peaks/seacr/{sample_id}.bed",
            seacr=temp("seqnado_output/peaks/seacr/{sample_id}_seacr.txt"),
            noM=temp("seqnado_output/bedgraphs/{sample_id}.nochrM.bedGraph"),
        log:
            "seqnado_output/logs/seacr/{sample_id}.log",
        params:
            threshold=CONFIG.third_party_tools.seacr.threshold,
            norm=CONFIG.third_party_tools.seacr.normalization,
            stringency=CONFIG.third_party_tools.seacr.stringency,
            prefix=lambda wc, output: Path(output.peaks).parent / Path(output.peaks).name,
        threads: 1
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=5, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        shell:
            """
            awk '$1 != "chrM"' {input.treatment} > {output.noM}
            SEACR_1.3.sh {output.noM} {params.threshold} {params.norm} {params.stringency} {output.peaks} > {log} 2>&1 || touch {params.prefix}.{params.stringency}.bed
            mv {params.prefix}.{params.stringency}.bed {output.seacr}
            cut -f 1-3 {output.seacr} > {output.peaks}
            """
        