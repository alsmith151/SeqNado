import os
import seqnado.utils


rule fastqc_raw:
    input:
        "seqnado_output/fastq/{sample}_{read}.fastq.gz",
    output:
        html="capcruncher_output/interim/qc/fastqc_raw/{sample}_{read}.html",
        zip="capcruncher_output/interim/qc/fastqc_raw/{sample}_{read}_fastqc.zip",  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra="--quiet",
    threads: 1
    resources:
        mem_mb=1500,
    log:
        "seqnado_output/logs/fastqc_raw/{sample}_{read}.log",
    wrapper:
        "v3.0.1/bio/fastqc"


use rule fastqc_raw as fastqc_trimmed with:
    input:
        "seqnado_output/trimmed/{sample}_{read}.fastq.gz",
    output:
        html="seqnado_output/qc/fastqc_trimmed/{sample}_{read}.html",
        zip="seqnado_output/qc/fastqc_trimmed/{sample}_{read}_fastqc.zip",  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        "seqnado_output/logs/fastqc_trimmed/{sample}_{read}.log",


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
            expand(
                "seqnado_output/qc/alignment_filtered/{sample}.txt",
                sample=SAMPLE_NAMES,
            ),
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
            expand(
                "seqnado_output/qc/alignment_filtered/{sample}.txt",
                sample=SAMPLE_NAMES,
            ),
        output:
            "seqnado_output/qc/full_qc_report.html",
        log:
            "seqnado_output/logs/multiqc.log",
        resources:
            mem_mb=lambda wildcards, attempt: 1000 * 2**attempt,
        shell:
            "multiqc -o seqnado_output/qc seqnado_output/qc -n full_qc_report.html --force > {log} 2>&1"
