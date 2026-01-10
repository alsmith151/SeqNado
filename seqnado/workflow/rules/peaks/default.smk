from seqnado import FileType, DataScalingTechnique, Assay
from seqnado.inputs import FastqCollection, FastqCollectionForIP
from seqnado.config.third_party_tools import CommandLineArguments


def get_control_file(wildcards, file_type: FileType):
    """
    Get the control file for a given treatment sample.

    Args:
        wildcards: Snakemake wildcards object containing 'sample' and 'treatment'.
        file_type (FileType): The type of file to retrieve (e.g., FileType.BAM, FileType.BIGWIG, FileType.TAG_DIRECTORY).

    Returns:
        str or None: The path to the control file if it exists, otherwise None.
    """

    if isinstance(INPUT_FILES, FastqCollection):
        return "NO_CONTROL_PRESENT_FOR_FILE"  # Dummy value to prevent the rule from ever being able to run if not an IP
    elif isinstance(INPUT_FILES, FastqCollectionForIP):
        search_term = wildcards.sample_id
        control_name = INPUT_FILES.get_control_performed(search_term)
        if not control_name:
            return "NO_CONTROL_PRESENT_FOR_FILE"

    # Match control with correct filetype
    match file_type:
        case FileType.BAM:
            return OUTPUT_DIR + f"/aligned/{control_name}.bam"
        case FileType.BIGWIG:
            return (
                OUTPUT_DIR
                + f"/bigwigs/deeptools/{DataScalingTechnique.UNSCALED.value}/{control_name}.bigWig"
            )
        case FileType.TAG_DIRECTORY:
            return OUTPUT_DIR + f"/tag_dirs/{control_name}"
        case _:
            raise ValueError(
                f"Unsupported file type '{file_type}' for control file retrieval"
            )


def correct_macs_options(wildcards, options: CommandLineArguments):

    # Correct the options based on whether the input is paired or not
    if hasattr(wildcards, "group"):
        sample = (
            SAMPLE_GROUPINGS.get_grouping("consensus")
            .get_group(wildcards.group)
            .samples
        )
        if not all(INPUT_FILES.is_paired_end(sample_name) for sample_name in sample):
            options = CommandLineArguments(value=options, exclude={"-f"})
    else:

        sample_id = wildcards.sample_id
        is_paired = INPUT_FILES.is_paired_end(sample_id)
        if not is_paired:
            options = CommandLineArguments(value=options.value, exclude={"-f"})
    return options


def get_lanceotron_call_peaks_threshold(wildcards):
    options = str(CONFIG.third_party_tools.lanceotron.call_peaks.command_line_arguments)
    threshold_pattern = re.compile(r"\-c\s+(\d+\.?\d*)")
    match = threshold_pattern.search(options)
    if match:
        return match.group(1)
    else:
        return "0.5"  # Default threshold


rule macs2_with_input:
    input:
        treatment=OUTPUT_DIR + "/aligned/{sample_id}.bam",
        control=lambda wc: get_control_file(wc, file_type=FileType.BAM),
    output:
        peaks=OUTPUT_DIR + "/peaks/macs2/{sample_id}.bed",
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=lambda wc: str(
            correct_macs_options(
                wc, CONFIG.third_party_tools.macs.call_peaks.command_line_arguments
            )
        ),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
        narrow_peak=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=6, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "docker://quay.io/biocontainers/macs2:2.1.1.20160309--py27r3.3.1_1"
    log:
        OUTPUT_DIR + "/logs/macs2/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/peaks/macs2/{sample_id}.tsv"
    message:
        "Calling peaks with MACS2 for sample {wildcards.sample_id}"
    shell:
        """
    if ! macs2 callpeak -t {input.treatment} -c {input.control} -n {params.basename} {params.options} > {log} 2>&1; then
        touch {output.peaks}
    else
        awk 'BEGIN{{OFS="\\t"}} !/^#/ && !/^chr[[:space:]]+start[[:space:]]+end/ && !/^$/ {{print $1, $2, $3}}' {params.narrow_peak} > {output.peaks} 2>> {log}
    fi
    """


rule macs2_no_input:
    input:
        treatment=OUTPUT_DIR + "/aligned/{sample_id}.bam",
    output:
        peaks=OUTPUT_DIR + "/peaks/macs2/{sample_id}.bed",
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=lambda wc: str(
            correct_macs_options(
                wc, CONFIG.third_party_tools.macs.call_peaks.command_line_arguments
            )
        ),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
        narrow_peak=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "docker://quay.io/biocontainers/macs2:2.1.1.20160309--py27r3.3.1_1"
    log:
        OUTPUT_DIR + "/logs/macs2/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/peaks/macs2/{sample_id}.tsv"
    message:
        "Calling peaks with MACS2 for sample {wildcards.sample_id}"
    shell:
        """
    if ! macs2 callpeak -t {input.treatment} -n {params.basename} {params.options} > {log} 2>&1; then
        touch {output.peaks}
    else
        awk 'BEGIN{{OFS="\\t"}} !/^#/ && !/^chr[[:space:]]+start[[:space:]]+end/ && !/^$/ {{print $1, $2, $3}}' {params.narrow_peak} > {output.peaks} 2>> {log}
    fi
    """


ruleorder: macs2_with_input > macs2_no_input


# MACS3 rules (same as MACS2 but with macs3 container)
rule macs3_with_input:
    input:
        treatment=OUTPUT_DIR + "/aligned/{sample_id}.bam",
        control=lambda wc: get_control_file(wc, file_type=FileType.BAM),
    output:
        peaks=OUTPUT_DIR + "/peaks/macs3/{sample_id}.bed",
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=lambda wc: str(
            correct_macs_options(
                wc, CONFIG.third_party_tools.macs.call_peaks.command_line_arguments
            )
        ),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
        narrow_peak=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=6, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "docker://quay.io/biocontainers/macs3:3.0.3--py39h0699b22_0"
    log:
        OUTPUT_DIR + "/logs/macs3/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/peaks/macs3/{sample_id}.tsv"
    message:
        "Calling peaks with MACS3 for sample {wildcards.sample_id}"
    shell:
        """
    if ! macs3 callpeak -t {input.treatment} -c {input.control} -n {params.basename} {params.options} > {log} 2>&1; then
        touch {output.peaks}
    else
        awk 'BEGIN{{OFS="\\t"}} !/^#/ && !/^chr[[:space:]]+start[[:space:]]+end/ && !/^$/ {{print $1, $2, $3}}' {params.narrow_peak} > {output.peaks} 2>> {log}
    fi
    """


rule macs3_no_input:
    input:
        treatment=OUTPUT_DIR + "/aligned/{sample_id}.bam",
    output:
        peaks=OUTPUT_DIR + "/peaks/macs3/{sample_id}.bed",
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=lambda wc: str(
            correct_macs_options(
                wc, CONFIG.third_party_tools.macs.call_peaks.command_line_arguments
            )
        ),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
        narrow_peak=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "docker://quay.io/biocontainers/macs3:3.0.3--py39h0699b22_0"
    log:
        OUTPUT_DIR + "/logs/macs3/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/peaks/macs3/{sample_id}.tsv"
    message:
        "Calling peaks with MACS3 for sample {wildcards.sample_id}"
    shell:
        """
    if ! macs3 callpeak -t {input.treatment} -n {params.basename} {params.options} > {log} 2>&1; then
        touch {output.peaks}
    else
        awk 'BEGIN{{OFS="\\t"}} !/^#/ && !/^chr[[:space:]]+start[[:space:]]+end/ && !/^$/ {{print $1, $2, $3}}' {params.narrow_peak} > {output.peaks} 2>> {log}
    fi
    """


ruleorder: macs3_with_input > macs3_no_input


rule homer_with_input:
    input:
        treatment=OUTPUT_DIR + "/tag_dirs/{sample_id}",
        control=lambda wc: get_control_file(wc, file_type=FileType.TAG_DIRECTORY),
    output:
        peaks=OUTPUT_DIR + "/peaks/homer/{sample_id}.bed",
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=str(CONFIG.third_party_tools.homer.find_peaks.command_line_arguments),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/homer/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/homer/{sample_id}.tsv"
    message:
        "Calling peaks with HOMER for sample {wildcards.sample_id}"
    shell:
        """
    findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp  -i {input.control} > {log} 2>&1 &&
    pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
    rm {output.peaks}.tmp
    """


rule homer_no_input:
    input:
        treatment=OUTPUT_DIR + "/tag_dirs/{sample_id}",
    output:
        peaks=OUTPUT_DIR + "/peaks/homer/{sample_id}.bed",
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=str(CONFIG.third_party_tools.homer.find_peaks.command_line_arguments),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/homer/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/homer/{sample_id}.tsv"
    message:
        "Calling peaks with HOMER for sample {wildcards.sample_id}"
    shell:
        """
    findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp > {log} 2>&1 &&
    pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
    rm {output.peaks}.tmp
    """


ruleorder: homer_with_input > homer_no_input


rule lanceotron_with_input:
    input:
        treatment=OUTPUT_DIR + "/bigwigs/deeptools/unscaled/{sample_id}.bigWig",
        control=lambda wc: get_control_file(wc, file_type=FileType.BIGWIG),
    output:
        peaks=OUTPUT_DIR + "/peaks/lanceotron/{sample_id}.bed",
        ltron_peaks=temp(OUTPUT_DIR + "/peaks/lanceotron/{sample_id}_L-tron.bed"),
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        threshold=lambda wildcards: get_lanceotron_call_peaks_threshold(wildcards),
        options=str(
            CONFIG.third_party_tools.lanceotron.call_peaks.command_line_arguments
        ),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    container:
        "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest"
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=12, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=6, attempts=attempt, scale=SCALE_RESOURCES
        ),
    log:
        OUTPUT_DIR + "/logs/lanceotron/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/lanceotron/{sample_id}.tsv"
    message:
        "Calling peaks with LanceOtron for sample {wildcards.sample_id}"
    shell:
        """
    lanceotron callPeaksInput {input.treatment} -i {input.control} -f {params.outdir} --skipheader > {log} 2>&1 &&
    cat {output.ltron_peaks} | awk 'BEGIN{{OFS="\\t"}} $4 >= {params.threshold} {{print $1, $2, $3}}' > {output.peaks} 
    """


rule lanceotron_no_input:
    input:
        treatment=OUTPUT_DIR + "/bigwigs/deeptools/unscaled/{sample_id}.bigWig",
    output:
        peaks=OUTPUT_DIR + "/peaks/lanceotron/{sample_id}.bed",
        ltron_peaks=temp(OUTPUT_DIR + "/peaks/lanceotron/{sample_id}_L-tron.bed"),
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=str(
            CONFIG.third_party_tools.lanceotron.call_peaks.command_line_arguments
        ),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    container:
        "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=12, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=6, attempts=attempt, scale=SCALE_RESOURCES
        ),
    log:
        OUTPUT_DIR + "/logs/lanceotron/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/lanceotron/{sample_id}.tsv"
    message:
        "Calling peaks with LanceOtron for sample {wildcards.sample_id}"
    shell:
        """
    lanceotron callPeaks {input.treatment} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
    cat {output.ltron_peaks} | cut -f 1-3 > {output.peaks}
    """


ruleorder: lanceotron_with_input > lanceotron_no_input



rule seacr:
    input:
        treatment=OUTPUT_DIR + "/bedgraphs/{sample_id}.bedGraph",
    output:
        peaks=OUTPUT_DIR + "/peaks/seacr/{sample_id}.bed",
        seacr=temp(OUTPUT_DIR + "/peaks/seacr/{sample_id}_seacr.txt"),
        noM=temp(OUTPUT_DIR + "/bedgraphs/{sample_id}.nochrM.bedGraph"),
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        threshold=CONFIG.third_party_tools.seacr.threshold,
        norm=CONFIG.third_party_tools.seacr.normalization,
        stringency=CONFIG.third_party_tools.seacr.stringency,
        prefix=lambda wc, output: Path(output.peaks).parent
        / Path(output.peaks).name,
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=5, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/seacr/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/seacr/{sample_id}.tsv"
    message:
        "Calling peaks with SEACR for sample {wildcards.sample_id}"
    shell:
        """
    awk '$1 != "chrM"' {input.treatment} > {output.noM}
    SEACR_1.3.sh {output.noM} {params.threshold} {params.norm} {params.stringency} {output.peaks} > {log} 2>&1 || touch {params.prefix}.{params.stringency}.bed
    mv {params.prefix}.{params.stringency}.bed {output.seacr}
    cut -f 1-3 {output.seacr} > {output.peaks}
    """
