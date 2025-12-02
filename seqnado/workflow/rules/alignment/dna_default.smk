from seqnado.helpers import define_time_requested, define_memory_requested
from seqnado import SpikeInMethod

rule align_paired:
    input:
        fq1=OUTPUT_DIR + "/trimmed/{sample}_1.fastq.gz",
        fq2=OUTPUT_DIR + "/trimmed/{sample}_2.fastq.gz",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/raw/{sample}.bam"),
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.align.command_line_arguments),
        rg="--rg-id {sample} --rg SM:{sample}",
    threads: CONFIG.third_party_tools.bowtie2.align.threads,
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/align/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}.tsv",
    message: "Aligning reads for sample {wildcards.sample} using Bowtie2",
    shell: """
    bowtie2 \
        -p {threads} \
        -x {params.index} \
        -1 {input.fq1} \
        -2 {input.fq2} \
        {params.rg} \
        {params.options} \
        2> {log} \
    | samtools view -bS - > {output.bam}
    """

rule align_single:
    input:
        fq1=OUTPUT_DIR + "/trimmed/{sample}.fastq.gz",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/raw/{sample}.bam"),
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.align.command_line_arguments),
        rg="--rg-id {sample} --rg SM:{sample}",
    threads: CONFIG.third_party_tools.bowtie2.align.threads,
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/align/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}.tsv",
    message: "Aligning reads for sample {wildcards.sample} using Bowtie2",
    shell: """
    bowtie2 \
        -p {threads} \
        -x {params.index} \
        -U {input.fq1} \
        {params.rg} \
        {params.options} \
        2> {log} \
    | samtools view -bS - > {output.bam}
    """


ruleorder: align_paired > align_single