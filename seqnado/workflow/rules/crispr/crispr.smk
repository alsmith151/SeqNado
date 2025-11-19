from seqnado.helpers import check_options, define_time_requested, define_memory_requested

rule crispr_trimming:
    input:
        fq=OUTPUT_DIR + "/fastqs/{sample}.fastq.gz",

    output:
        trimmed=temp(OUTPUT_DIR + "/trimmed/{sample}.fastq.gz"),
    params:
        options=check_options(config["cutadapt"]["options"]),
        trim_dir=OUTPUT_DIR + "/trimmed",
    threads: config["cutadapt"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/trimming/{sample}.log",
    shell:"""
    cutadapt {params.options} -o {output.trimmed} {input.fq} > {log} 2>&1
    """


rule align_crispr:
    input:
        fq1=OUTPUT_DIR + "/trimmed/{sample}.fastq.gz",
    params:
        index=config["genome"]["index"],
        options=check_options(config["bowtie2"]["options"]),
    output:
        bam=temp(OUTPUT_DIR + "/aligned/raw/{sample}.bam"),
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: config["bowtie2"]["threads"]
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/align/{sample}.log",
    shell:"""
    bowtie2 -p {threads} -x {params.index} -U {input.fq1} {params.options} 2> {log} |
    samtools view -bS - > {output.bam}
    """

rule feature_counts:
    input:
        bam=expand(OUTPUT_DIR + "/aligned/{sample}.bam", sample=SAMPLE_NAMES),
        bai=expand(OUTPUT_DIR + "/aligned/{sample}.bam.bai", sample=SAMPLE_NAMES),
        annotation=config["genome"]["gtf"],
    output:
        counts=OUTPUT_DIR + "/readcounts/feature_counts/read_counts.tsv",
    params:
        options=check_options(config["featurecounts"]["options"]),
    threads: config["featurecounts"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/readcounts/featurecounts/featurecounts.log",
    shell:"""
    featureCounts \
    -a \
    {input.annotation} \
    -T \
    {threads} \
    {params.options} \
    -o \
    {output.counts} \
    {input.bam} \
    > {log} 2>&1
    """