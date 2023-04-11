from typing import Literal
import seqnado.utils
  
rule macs2_with_input:
    input:
        unpack(lambda wc: seqnado.utils.pair_treatment_and_control_for_peak_calling(wc, samples=FASTQ_SAMPLES, assay=ASSAY, filetype="bam")),
    output:
        peaks="seqnado_output/peaks/macs/{sample}_{antibody}.bed",
    params:
        options=seqnado.utils.check_options(config["macs"]["callpeak"]),
        narrow=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        "seqnado_output/logs/macs/{sample}_{antibody}.log",
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} -n seqnado_output/peaks/macs/{wildcards.treatment} -f BAMPE {params.options} > {log} 2>&1 &&
        cat {params.narrow} | cut -f 1-3 > {output.peaks}
        """

rule macs2_no_input:
    input:
        treatment="seqnado_output/aligned/{sample}_{antibody}.bam",
    output:
        peaks="seqnado_output/peaks/macs/{sample}_{antibody}.bed",
    params:
        options=seqnado.utils.check_options(config["macs"]["callpeak"]),
        narrow=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        "seqnado_output/logs/macs/{sample}_{antibody}.log",
    shell:
        """
        macs2 callpeak -t {input.treatment} -n {params.basename} -f BAMPE {params.options} > {log} 2>&1 &&
        cat {params.narrow} | cut -f 1-3 > {output.peaks}
        """


rule homer_with_input:
    input:
        unpack(lambda wc: seqnado.utils.pair_treatment_and_control_for_peak_calling(wc, samples=FASTQ_SAMPLES, assay=ASSAY, filetype="tag")),
    output:
        peaks="seqnado_output/peaks/homer/{sample}_{antibody}.bed",
    log:
        "seqnado_output/logs/homer/findPeaks/{sample}_{antibody}.log",
    params:
        options=seqnado.utils.check_options(config["homer"]["findpeaks"]),
    threads: 1
    resources:
        mem_mb=2000,
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp  -i {input.control} > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule homer_no_input:
    input:
        treatment="seqnado_output/tag_dirs/{sample}_{antibody}/",
    output:
        peaks="seqnado_output/peaks/homer/{sample}_{antibody}.bed",
    log:
        "seqnado_output/logs/homer/findPeaks_{sample}_{antibody}.log",
    params:
        options=seqnado.utils.check_options(config["homer"]["findpeaks"]),
    threads: 1
    resources:
        mem_mb=1024,
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule lanceotron_with_input:
    input:
        unpack(lambda wc: seqnado.utils.pair_treatment_and_control_for_peak_calling(wc, samples=FASTQ_SAMPLES, assay=ASSAY, filetype="bigwig")),
    output:
        peaks="seqnado_output/peaks/lanceotron/{sample}_{antibody}.bed",
    log:
        "seqnado_output/logs/lanceotron/{sample}_{antibody}.log",
    params:
        options=seqnado.utils.check_options(config["lanceotron"]["callpeak"]),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
    container:
        "library://asmith151/seqnado/seqnado_extra:v1"
    threads: 1
    resources:
        mem_mb=1024 * 10,
    shell:
        """
        lanceotron callPeaksInput {input.treatment} -i {input.control} -f {params.outdir} {params.options} --skipheader > {log} 2>&1 &&
        cat {params.outdir}/{wildcards.sample}_{wildcards.antibody}_L-tron.bed | cut -f 1-3 > {output.peaks}
        """


rule lanceotron_no_input:
    input:
        treatment="seqnado_output/bigwigs/deeptools/{sample}_{antibody}.bigWig",
    output:
        peaks="seqnado_output/peaks/lanceotron/{sample}_{antibody}.bed",
    log:
        "seqnado_output/logs/lanceotron/{sample}_{antibody}.log",
    params:
        options=seqnado.utils.check_options(config["lanceotron"]["callpeak"]),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
    threads: 1
    container:
        "library://asmith151/seqnado/seqnado_extra:v1"
    resources:
        mem_mb=1024 * 10,
    shell:
        """
        lanceotron callPeaks {input.treatment} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
        cat {params.outdir}/{wildcards.sample}_{wildcards.antibody}_L-tron.bed | cut -f 1-3 > {output.peaks}
        """


ruleorder: lanceotron_with_input > lanceotron_no_input
ruleorder: homer_with_input > homer_no_input
ruleorder: macs2_with_input > macs2_no_input
