from seqnado.helpers import define_time_requested, define_memory_requested



rule align_paired:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.command_line_arguments),
        rg="--rg-id {sample} --rg SM:{sample}",
    output:
        bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
    threads: config["bowtie2"]["threads"]
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/align/{sample}.log",
    shell:
        """
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
        fq1="seqnado_output/trimmed/{sample}.fastq.gz",
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.command_line_arguments),
        rg="--rg-id {sample} --rg SM:{sample}",
    threads: config["bowtie2"]["threads"],
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/align/{sample}.log",
    output:
        bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
    shell:
        """
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
