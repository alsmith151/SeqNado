


rule split_reads_aligned_to_viewpoints:
    input:
        bam=OUTPUT_DIR + "/aligned/aligned_to_viewpoints/{sample}.bam",
    output:
        fq=temp(OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.sliced.fastq.gz"),
    threads: 1
    resources:
        mem="1GB",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/split_reads/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/split_reads/{sample}.tsv",
    message: "Splitting reads aligned to viewpoints for sample {wildcards.sample}",
    script:
        "../../scripts/mcc_split_reads_aligned_to_viewpoints.py"


rule align_mcc_reads_to_genome:
    input:
        fq1=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.sliced.fastq.gz",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/initial_alignment/{sample}.bam"),
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.align.command_line_arguments),
        rg="--rg-id {sample} --rg SM:{sample}",
    threads: CONFIG.third_party_tools.bowtie2.align.threads,
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/align_mcc/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align_mcc/{sample}.tsv",
    message: "Aligning MCC reads to genome for sample {wildcards.sample}",
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


rule align_unmapped_reads_to_genome:
    input:
        bam=OUTPUT_DIR + "/aligned/initial_alignment/{sample}.bam",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/second_alignment/{sample}.bam"),
        bai=temp(OUTPUT_DIR + "/aligned/second_alignment/{sample}.bam.bai"),
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.align.command_line_arguments),
    threads: CONFIG.third_party_tools.bowtie2.align.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/realign/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/realign/{sample}.tsv",
    message: "Aligning unmapped MCC reads to genome for sample {wildcards.sample}",
    shell: """
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
        bam=temp(OUTPUT_DIR + "/aligned/raw/{sample}.bam"),
    threads: CONFIG.third_party_tools.bowtie2.align.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/combine/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/combine/{sample}.tsv",
    message: "Combining genome-mapped reads for sample {wildcards.sample}",
    shell: """
    samtools merge -@ {threads} {output.bam} {input.bam1} {input.bam2} &&
    samtools view -F 4 -b {output.bam} > {output.bam}.tmp &&
    mv {output.bam}.tmp {output.bam} &&
    samtools sort -n -@ {threads} -o {output.bam}_sorted {output.bam} &&
    mv {output.bam}_sorted {output.bam}
    """