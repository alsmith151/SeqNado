from typing import Literal
import seqnado.utils
  

rule macs2_no_input:
    input:
        treatment="seqnado_output/aligned/{sample}.bam",
    output:
        peaks="seqnado_output/peaks/macs/{sample}.bed",
    params:
        options=seqnado.utils.check_options(config["macs"]["callpeak"]),
        narrow=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        "seqnado_output/logs/macs/{sample}.log",
    shell:
        """
        macs2 callpeak -t {input.treatment} -n {params.basename} -f BAMPE {params.options} > {log} 2>&1 &&
        cat {params.narrow} | cut -f 1-3 > {output.peaks}
        """

rule homer_no_input:
    input:
        treatment="seqnado_output/tag_dirs/{sample}/",
    output:
        peaks="seqnado_output/peaks/homer/{sample}.bed",
    log:
        "seqnado_output/logs/homer/findPeaks_{sample}.log",
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

rule lanceotron_no_input:
    input:
        treatment="seqnado_output/bigwigs/deeptools/{sample}.bigWig",
    output:
        peaks="seqnado_output/peaks/lanceotron/{sample}.bed",
    log:
        "seqnado_output/logs/lanceotron/{sample}.log",
    params:
        options=seqnado.utils.check_options(config["lanceotron"]["callpeak"]),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
    threads: 1
    resources:
        mem_mb=1024 * 10,
    shell:
        """
        lanceotron callPeaks {input.treatment} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
        cat {params.outdir}/{wildcards.sample}_L-tron.bed | cut -f 1-3 > {output.peaks}
        """
