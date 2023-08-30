import os
import seqnado.utils

rule fastqc_raw:
    input:
        unpack(lambda wc: seqnado.utils.translate_fq_files(wc, samples=FASTQ_SAMPLES, paired=False)),
    output:
        qc="seqnado_output/qc/fastqc_raw/{sample}_{read}_fastqc.html",
    params:
        outdir="seqnado_output/qc/fastqc_raw",
        tmpdir="seqnado_output/qc/fastqc_raw/{sample}_{read}",
        basename=lambda wc, output: seqnado.utils.get_fq_filestem(wc, samples=FASTQ_SAMPLES),
    threads: 1
    resources:
        mem_mb=1500,
    log:
        "seqnado_output/logs/fastqc_raw/{sample}_{read}.log",
    shell:
        """
        mkdir -p {params.tmpdir} &&
        fastqc -o {params.tmpdir} {input.fq} > {log} 2>&1 &&
        mv {params.tmpdir}/{params.basename}_fastqc.html {output.qc} &&
        mv {params.tmpdir}/{params.basename}_fastqc.zip {params.outdir}/{wildcards.sample}_{wildcards.read}_fastqc.zip &&
        rm -r {params.tmpdir}
        """

rule fastqc_trimmed:
    input:
        fq="seqnado_output/trimmed/{sample}_{read}.fastq.gz",
    output:
        qc="seqnado_output/qc/fastqc_trimmed/{sample}_{read}_fastqc.html",
    params:
        outdir="seqnado_output/qc/fastqc_trimmed",
    log:
        "seqnado_output/logs/fastqc_trimmed/{sample}_{read}.log",
    resources:
        mem_mb=1500,
    threads: 4
    shell:
        """fastqc -o {params.outdir} {input.fq} > {log} 2>&1"""


rule samtools_stats:
    input:
        bam="seqnado_output/aligned/raw/{sample}.bam",
    output:
        stats="seqnado_output/qc/alignment_raw/{sample}.txt",
    threads: 1
    resources:
        mem_mb=1000,
    shell:
        """samtools stats {input.bam} > {output.stats}"""


use rule samtools_stats as samtools_stats_filtered with:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
    output:
        stats="seqnado_output/qc/alignment_filtered/{sample}.txt",

if config["split_fastq"] == "False":
    rule multiqc:
        input:
            expand(
                "seqnado_output/qc/fastqc_raw/{sample}_{read}_fastqc.html",
                sample=SAMPLE_NAMES,
                read=[1, 2],
            ),
            expand(
                "seqnado_output/qc/fastqc_trimmed/{sample}_{read}_fastqc.html",
                sample=SAMPLE_NAMES,
                read=[1, 2],
            ),
            expand("seqnado_output/qc/alignment_raw/{sample}.txt", sample=SAMPLE_NAMES),
            expand("seqnado_output/qc/alignment_filtered/{sample}.txt", sample=SAMPLE_NAMES),
        output:
            "seqnado_output/qc/full_qc_report.html",
        log:
            "seqnado_output/logs/multiqc.log",
        resources:
            mem_mb=lambda wildcards, attempt: 2000 * 2**attempt,
        shell:
            "multiqc -o seqnado_output/qc seqnado_output/qc -n full_qc_report.html --force > {log} 2>&1"
else:
    rule multiqc:
        input:
            expand(
                "seqnado_output/qc/fastqc_raw/{sample}_{read}_fastqc.html",
                sample=SAMPLE_NAMES,
                read=[1, 2],
            ),
            expand("seqnado_output/qc/alignment_raw/{sample}.txt", sample=SAMPLE_NAMES),
            expand("seqnado_output/qc/alignment_filtered/{sample}.txt", sample=SAMPLE_NAMES),
        output:
            "seqnado_output/qc/full_qc_report.html",
        log:
            "seqnado_output/logs/multiqc.log",
        resources:
            mem_mb=lambda wildcards, attempt: 1000 * 2**attempt,
        shell:
            "multiqc -o seqnado_output/qc seqnado_output/qc -n full_qc_report.html --force > {log} 2>&1"

    