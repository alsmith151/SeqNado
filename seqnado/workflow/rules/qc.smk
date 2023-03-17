def get_fq_file(wc):
    return {"fq": FASTQ_SAMPLES.translation[f"{wc.sample}_{wc.read}.fastq.gz"]}

rule fastqc_raw:
    input:
        unpack(get_fq_file),
    output:
        qc="seqnado_output/qc/fastq_raw/{sample}_{read}_fastqc.html",
    params:
        outdir="seqnado_output/qc/fastq_raw",
    threads: 4
    log:
        "seqnado_output/logs/fastqc_raw/{sample}_{read}.log",
    shell:
        "fastqc -o {params.outdir} {input.fq} > {log} 2>&1"


use rule fastqc_raw as fastqc_trimmed with:
    input:
        fq="seqnado_output/trimmed/{sample}_{read}.fastq.gz",
    output:
        qc="seqnado_output/qc/fastq_trimmed/{sample}_{read}_fastqc.html",
    params:
        outdir="seqnado_output/qc/fastq_trimmed",
    log:
        "seqnado_output/logs/fastqc_trimmed/{sample}_{read}.log",


rule samtools_stats:
    input:
        bam="seqnado_output/aligned/sorted/{sample}.bam",
    output:
        stats="seqnado_output/qc/alignment_raw/{sample}.txt",
    threads: 1
    shell:
        """samtools stats {input.bam} > {output.stats}"""


use rule samtools_stats as samtools_stats_filtered with:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
    output:
        stats="seqnado_output/qc/alignment_filtered/{sample}.txt",


rule multiqc:
    input:
        expand(
            "seqnado_output/qc/fastq_raw/{sample}_{read}_fastqc.html",
            sample=SAMPLE_NAMES,
            read=[1, 2],
        ),
        expand(
            "seqnado_output/qc/fastq_trimmed/{sample}_{read}_fastqc.html",
            sample=SAMPLE_NAMES,
            read=[1, 2],
        ),
        expand("seqnado_output/qc/alignment_raw/{sample}.txt", sample=SAMPLE_NAMES),
        expand("seqnado_output/qc/alignment_filtered/{sample}.txt", sample=SAMPLE_NAMES),
    output:
        "seqnado_output/qc/full_qc_report.html",
    log:
        "seqnado_output/logs/multiqc.log",
    shell:
        "multiqc -o qc qc -n full_qc_report.html --force > {log} 2>&1"
