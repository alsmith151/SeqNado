from typing import Literal
from seqnado.helpers import check_options, define_time_requested, define_memory_requested
import re


def get_lanceotron_threshold(wildcards):
    options = config["lanceotron"]["callpeak"]
    threshold_pattern = re.compile(r"\-c\s+(\d+.?\d*)")
    threshold = threshold_pattern.search(options).group(1)
    return threshold


rule macs2_no_input:
    input:
        treatment="seqnado_output/aligned/{sample}.bam",
    output:
        peaks="seqnado_output/peaks/macs/{sample}.bed",
    params:
        options=check_options(config["macs"]["callpeak"]),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/macs/{sample}.log",
    container:
        "docker://quay.io/biocontainers/macs2:2.2.9.1--py311haab0aaa_3"
    shell:
        """
        macs2 callpeak -t {input.treatment} -n {params.basename} -f BAMPE {params.options} > {log} 2>&1 &&
        cat {params.raw} | grep -v '^#' | grep -vE '^chr\\s+start\\s+end.*' | grep -v '^$' | cut -f 1-3 > {output.peaks}
        """


rule homer_no_input:
    input:
        treatment="seqnado_output/tag_dirs/{sample}",
    output:
        peaks="seqnado_output/peaks/homer/{sample}.bed",
    log:
        "seqnado_output/logs/homer/{sample}.log",
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


rule lanceotron_no_input:
    input:
        treatment="seqnado_output/bigwigs/deeptools/unscaled/{sample}.bigWig",
    output:
        peaks="seqnado_output/peaks/lanceotron/{sample}.bed",
        ltron_peaks=temp("seqnado_output/peaks/lanceotron/{sample}_L-tron.bed"),
    log:
        "seqnado_output/logs/lanceotron/{sample}.log",
    params:
        options=check_options(config["lanceotron"]["callpeak"]),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
    threads: 1
    container:
        "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest"
    resources:
        mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt)}GB",
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    shell:"""
    lanceotron callPeaks {input.treatment} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
    cat {output.ltron_peaks} | cut -f 1-3 > {output.peaks}
    """
