from typing import Literal
import seqnado.utils
import re


def get_lanceotron_threshold(wildcards):
    options = config["lanceotron"]["callpeak"]
    threshold_pattern = re.compile(r"\-c\s+(\d+.?\d*)")
    threshold = threshold_pattern.search(options).group(1)
    return threshold

    
rule macs2_with_input:
    input:
        treatment = lambda wc: seqnado.utils.get_treatment_file(wc,assay=ASSAY, filetype="bam"),
        control = lambda wc: seqnado.utils.get_control_file(wc, design=DESIGN, assay=ASSAY, filetype="bam"),
    output:
        peaks="seqnado_output/peaks/macs/{wildcards.treatment}.bed",
    params:
        options=seqnado.utils.check_options(config["macs"]["callpeak"]),
        narrow=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
    threads: 1
    resources:
        mem_mb=2000,
        time='02:00:00',
    log:
        "seqnado_output/logs/macs/{wildcards.treatment}.bed",
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} -n seqnado_output/peaks/macs/{wildcards.treatment} -f BAMPE {params.options} > {log} 2>&1 &&
        cat {params.narrow} | cut -f 1-3 > {output.peaks}
        """

rule macs2_no_input:
    input:
        treatment = lambda wc: seqnado.utils.get_treatment_file(wc, assay=ASSAY, filetype="bam"),
    output:
        peaks="seqnado_output/peaks/macs/{treatment}.bed",
    params:
        options=seqnado.utils.check_options(config["macs"]["callpeak"]),
        narrow=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem_mb=2000,
        time='02:00:00',
    log:
        "seqnado_output/logs/macs/{treatment}.bed",
    shell:
        """
        macs2 callpeak -t {input.treatment} -n {params.basename} -f BAMPE {params.options} > {log} 2>&1 &&
        cat {params.narrow} | cut -f 1-3 > {output.peaks}
        """


rule homer_with_input:
    input:
        treatment=lambda wc: seqnado.utils.get_treatment_file(wc,assay=ASSAY, filetype="tag"),
        control=lambda wc: seqnado.utils.get_control_file(wc, design=DESIGN, assay=ASSAY, filetype="tag"),
    output:
        peaks="seqnado_output/peaks/homer/{treatment}.bed",
    log:
        "seqnado_output/logs/homer/{treatment}.bed",
    params:
        options=seqnado.utils.check_options(config["homer"]["findpeaks"]),
    threads: 1
    resources:
        mem_mb=2000,
        time='02:00:00',
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp  -i {input.control} > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule homer_no_input:
    input:
        treatment=lambda wc: seqnado.utils.get_treatment_file(wc, assay=ASSAY, filetype="tag"),
    output:
        peaks="seqnado_output/peaks/homer/{treatment}.bed",
    log:
        "seqnado_output/logs/homer/{treatment}.bed",
    params:
        options=seqnado.utils.check_options(config["homer"]["findpeaks"]),
    threads: 1
    resources:
        mem_mb=1024,
        time='02:00:00',
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule lanceotron_with_input:
    input:
        treatment=lambda wc: seqnado.utils.get_treatment_file(wc,assay=ASSAY, filetype="bigwig"),
        control=lambda wc: seqnado.utils.get_control_file(wc, design=DESIGN, assay=ASSAY, filetype="bigwig"),
    output:
        peaks="seqnado_output/peaks/lanceotron/{treatment}.bed",
    log:
        "seqnado_output/logs/lanceotron/{treatment}.bed",
    params:
        threshold=get_lanceotron_threshold,
        outdir=lambda wc, output: os.path.dirname(output.peaks),
    container:
        "library://asmith151/seqnado/seqnado_extra:latest"
    threads: 1
    resources:
        mem_mb=1024 * 10,
        time='04:00:00',
    shell:
        """
        lanceotron callPeaksInput {input.treatment} -i {input.control} -f {params.outdir} --skipheader > {log} 2>&1 &&
        cat {params.outdir}/{wildcards.treatment}_L-tron.bed | awk 'BEGIN{{OFS="\\t"}} $4 >= {params.threshold} {{print $1, $2, $3}}' > {output.peaks}
        """


rule lanceotron_no_input:
    input:
        treatment=lambda wc: seqnado.utils.get_treatment_file(wc,assay=ASSAY, filetype="bigwig"),
    output:
        peaks="seqnado_output/peaks/lanceotron/{treatment}.bed",
    log:
        "seqnado_output/logs/lanceotron/{treatment}.bed",
    params:
        options=seqnado.utils.check_options(config["lanceotron"]["callpeak"]),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
    threads: 1
    container:
        "library://asmith151/seqnado/seqnado_extra:latest"
    resources:
        mem_mb=1024 * 10,
        time='04:00:00',
    shell:
        """
        lanceotron callPeaks {input.treatment} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
        cat {params.outdir}/{wildcards.treatment}_L-tron.bed | cut -f 1-3 > {output.peaks}
        """


ruleorder: lanceotron_with_input > lanceotron_no_input > homer_with_input > homer_no_input > macs2_with_input > macs2_no_input
