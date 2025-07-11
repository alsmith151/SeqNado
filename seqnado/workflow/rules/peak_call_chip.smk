from typing import Literal
from seqnado.helpers import check_options, define_time_requested, define_memory_requested
import re
import pathlib


def get_lanceotron_threshold(wildcards):
    options = config["lanceotron"]["callpeak"]
    threshold_pattern = re.compile(r"\-c\s+(\d+.?\d*)")
    threshold = threshold_pattern.search(options).group(1)
    return threshold

def format_macs_options(wildcards, options):
    query_name = f"{wildcards.sample}_{wildcards.treatment}"    
    is_paired = DESIGN.query(query_name).is_paired
    options = check_options(options)
    if not is_paired:
        options = re.sub(r"-f BAMPE", "", options)
    if not options:
        return ""
    else:
        return options


def get_control_file(wildcards, file_type: Literal["bam", "tag", "bigwig"], allow_null=False):
    control_info = DESIGN.query(sample_name=f"{wildcards.sample}_{wildcards.treatment}", full_experiment=True)["control"]
    if control_info:
        match file_type:
            case "bam":
                return f"seqnado_output/aligned/{control_info.sample_name}.bam"
            case "tag":
                return f"seqnado_output/tag_dirs/{control_info.sample_name}"
            case "bigwig":
                return f"seqnado_output/bigwigs/deeptools/unscaled/{control_info.sample_name}.bigWig"
    else:
        if allow_null:
            return []
        else:
            return "UNDEFINED"
        if allow_null:
            return []
        else:
            return "UNDEFINED"

rule macs2_with_input:
    input:
        treatment="seqnado_output/aligned/{sample}_{treatment}.bam",
        control=lambda wc: get_control_file(wc, file_type="bam", allow_null=False),
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
        "docker://quay.io/biocontainers/macs2:latest"
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} -n {params.basename} {params.options} > {log} 2>&1 &&
        cat {params.raw} | grep -v '^#' | grep -vE '^chr\\s+start\\s+end.*' | grep -v '^$' | cut -f 1-3 > {output.peaks}
        """


rule macs2_no_input:
    input:
        treatment="seqnado_output/aligned/{sample}_{treatment}.bam",
        control=lambda wc: get_control_file(wc, file_type="bam", allow_null=True), 
    output:
        peaks="seqnado_output/peaks/macs/{sample}_{treatment}.bed",
    params:
        options=lambda wc: format_macs_options(wc, config["macs"]["callpeak"]),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/macs/{sample}_{treatment}.log",
    container:
        "docker://quay.io/biocontainers/macs2:latest"
    shell:
        """
        macs2 callpeak -t {input.treatment} -n {params.basename} {params.options} > {log} 2>&1 &&
        cat {params.raw} | grep -v '^#' | grep -vE '^chr\\s+start\\s+end.*' | grep -v '^$' | cut -f 1-3 > {output.peaks}
        """


rule homer_with_input:
    input:
        treatment="seqnado_output/tag_dirs/{sample}_{treatment}",
        control=lambda wc: get_control_file(wc, file_type="tag", allow_null=False),
    output:
        peaks="seqnado_output/peaks/homer/{sample}_{treatment}.bed",
    log:
        "seqnado_output/logs/homer/{sample}_{treatment}.log",
    params:
        options=check_options(config["homer"]["findpeaks"]),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp  -i {input.control} > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule homer_no_input:
    input:
        treatment="seqnado_output/tag_dirs/{sample}_{treatment}",
        control=lambda wc: get_control_file(wc, file_type="tag", allow_null=True),
    output:
        peaks="seqnado_output/peaks/homer/{sample}_{treatment}.bed",
    log:
        "seqnado_output/logs/homer/{sample}_{treatment}.log",
    params:
        options=check_options(config["homer"]["findpeaks"]),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule lanceotron_with_input:
    input:
        treatment="seqnado_output/bigwigs/deeptools/unscaled/{sample}_{treatment}.bigWig",
        control=lambda wc: get_control_file(wc, file_type="bigwig", allow_null=False),
    output:
        peaks="seqnado_output/peaks/lanceotron/{sample}_{treatment}.bed",
        ltron_peaks=temp("seqnado_output/peaks/lanceotron/{sample}_{treatment}_L-tron.bed"),
    log:
        "seqnado_output/logs/lanceotron/{sample}_{treatment}.log",
    params:
        threshold=get_lanceotron_threshold,
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
        treatment="seqnado_output/bigwigs/deeptools/unscaled/{sample}_{treatment}.bigWig",
        control=lambda wc: get_control_file(wc, file_type="bigwig", allow_null=True),
    output:
        peaks="seqnado_output/peaks/lanceotron/{sample}_{treatment}.bed",
        ltron_peaks=temp("seqnado_output/peaks/lanceotron/{sample}_{treatment}_L-tron.bed"),
    log:
        "seqnado_output/logs/lanceotron/{sample}_{treatment}.log",
    params:
        options=check_options(config["lanceotron"]["callpeak"]),
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

rule seacr:
    input:
        treatment="seqnado_output/bedgraphs/{sample}_{treatment}.bedGraph",
    output:
        peaks="seqnado_output/peaks/seacr/{sample}_{treatment}.bed",
        seacr=temp("seqnado_output/peaks/seacr/{sample}_{treatment}_seacr.txt"),
        noM=temp("seqnado_output/bedgraphs/{sample}_{treatment}.nochrM.bedGraph"),
    log:
        "seqnado_output/logs/seacr/{sample}_{treatment}.log",
    params:
        threshold=config["seacr"].get("threshold", 0.01),
        norm=config["seacr"].get("norm", "non"),
        stringency=config["seacr"].get("stringency", "stringent"),
        prefix=lambda wc, output: pathlib.Path(output.peaks).parent / pathlib.Path(output.peaks).name,
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=5, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        """
        awk '$1 != "chrM"' {input.treatment} > {output.noM}
        SEACR_1.3.sh {output.noM} {params.threshold} {params.norm} {params.stringency} {output.peaks} > {log} 2>&1 || touch {params.prefix}.{params.stringency}.bed
        mv {params.prefix}.{params.stringency}.bed {output.seacr}
        cut -f 1-3 {output.seacr} > {output.peaks}
        """
    

ruleorder: lanceotron_with_input > lanceotron_no_input
ruleorder: homer_with_input > homer_no_input
ruleorder: macs2_with_input > macs2_no_input