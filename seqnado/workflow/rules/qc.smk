def get_fq_file(wc):
    return {"fq": FASTQ_SAMPLES.translation[f"{wc.sample}_{wc.read}.fastq.gz"]}

rule fastqc_raw:
    input:
        unpack(get_fq_file),
    output:
        qc="qc/fastq_raw/{sample}_{read}_fastqc.html",
    params:
        outdir="qc/fastq_raw",
    threads: 4
    log:
        "logs/fastqc_raw/{sample}_{read}.log",
    shell:
        "fastqc -o {params.outdir} {input.fq} > {log} 2>&1"


use rule fastqc_raw as fastqc_trimmed with:
    input:
        fq="trimmed/{sample}_{read}.fastq.gz",
    output:
        qc="qc/fastq_trimmed/{sample}_{read}_fastqc.html",
    params:
        outdir="qc/fastq_trimmed",
    log:
        "logs/fastqc_trimmed/{sample}_{read}.log",


rule samtools_stats:
    input:
        bam="aligned/sorted/{sample}.bam",
    output:
        stats="qc/alignment_raw/{sample}.txt",
    threads: 1
    shell:
        """samtools stats {input.bam} > {output.stats}"""


use rule samtools_stats as samtools_stats_filtered with:
    input:
        bam="aligned/{sample}.bam",
    output:
        stats="qc/alignment_filtered/{sample}.txt",


rule multiqc:
    input:
        expand(
            "qc/fastq_raw/{sample}_{read}_fastqc.html",
            sample=SAMPLE_NAMES,
            read=[1, 2],
        ),
        expand(
            "qc/fastq_trimmed/{sample}_{read}_fastqc.html",
            sample=SAMPLE_NAMES,
            read=[1, 2],
        ),
        expand("qc/alignment_raw/{sample}.txt", sample=SAMPLE_NAMES),
        expand("qc/alignment_filtered/{sample}.txt", sample=SAMPLE_NAMES),
    output:
        "qc/full_qc_report.html",
    log:
        "logs/multiqc.log",
    shell:
        "multiqc -o qc qc -n full_qc_report.html --force > {log} 2>&1"
