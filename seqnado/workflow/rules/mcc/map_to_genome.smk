
rule split_reads_aligned_to_viewpoints:
    input:
        bam="seqnado_output/aligned/aligned_to_viewpoints/{sample}.bam",
    output:
        fq="seqnado_output/mcc/{sample}/{sample}.sliced.fastq.gz",
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/split_reads/{sample}.log",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    script:
        "../scripts/mcc_split_reads_aligned_to_viewpoints.py"


use rule align_single as align_mcc_reads_to_genome with:
    input:
        fq1="seqnado_output/mcc/{sample}/{sample}.sliced.fastq.gz", 
    output:
        bam=temp("seqnado_output/aligned/initial_alignment/{sample}.bam"),


rule align_unmapped_reads_to_genome:
    input:
        bam="seqnado_output/aligned/initial_alignment/{sample}.bam",
    output:
        bam=temp("seqnado_output/aligned/second_alignment/{sample}.bam"),
        bai=temp("seqnado_output/aligned/second_alignment/{sample}.bam.bai"),
    threads: CONFIG.third_party_tools.bowtie2.align.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/realign/{sample}.log",
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.align.command_line_arguments),

    shell:
        """
        samtools view -b -f 4 {input.bam} | bowtie2 -p {threads} -x {params.index} -b - --very-sensitive-local 2>> {log} |
        samtools view -bS - > {output.bam} &&
        samtools sort -@ {threads} -o {output.bam}_sorted {output.bam} &&
        mv {output.bam}_sorted {output.bam} &&
        samtools index {output.bam}
        """

rule combine_genome_mapped_reads:
    input:
        bam1=rules.align_mcc_reads_to_genome.output.bam,
        bam2=rules.align_unmapped_reads_to_genome.output.bam,
    output:
        bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
    threads: CONFIG.third_party_tools.bowtie2.align.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/combine/{sample}.log",
    shell:
        """
        samtools merge -@ {threads} {output.bam} {input.bam1} {input.bam2} &&
        samtools view -F 4 -b {output.bam} > {output.bam}.tmp &&
        mv {output.bam}.tmp {output.bam} &&
        samtools sort -n -@ {threads} -o {output.bam}_sorted {output.bam} &&
        mv {output.bam}_sorted {output.bam}
        """