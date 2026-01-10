from seqnado.workflow.helpers.common import (
    define_memory_requested, 
    define_time_requested,
    get_alignment_input,
)


rule align_paired:
    input:
        fq1=lambda wildcards: get_alignment_input(wildcards, OUTPUT_DIR, CONFIG, paired=True)["fq1"],
        fq2=lambda wildcards: get_alignment_input(wildcards, OUTPUT_DIR, CONFIG, paired=True)["fq2"],
    params:
        index=str(CONFIG.genome.index.prefix),
        options=str(CONFIG.third_party_tools.star.align.command_line_arguments),
        prefix=OUTPUT_DIR + "/aligned/star/{sample}_",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/star/{sample}_Aligned.sortedByCoord.out.bam"),
        bam2=temp(
            OUTPUT_DIR + "/aligned/star/{sample}_Aligned.toTranscriptome.out.bam"
        ),
        log_out=temp(OUTPUT_DIR + "/aligned/star/{sample}_Log.final.out"),
    threads: CONFIG.third_party_tools.star.align.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=35, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/align/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}.tsv",
    message: "Aligning reads for sample {wildcards.sample} using STAR",
    shell: """
    STAR \
    --genomeDir {params.index} \
    --readFilesIn {input.fq1} {input.fq2} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN {threads} \
    --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} \
    --outFileNamePrefix {params.prefix} \
    {params.options} \
    > {log} 2>&1
    """


rule align_single:
    input:
        fq1=lambda wildcards: get_alignment_input(wildcards, OUTPUT_DIR, CONFIG, paired=False),
    params:
        index=str(CONFIG.genome.index.prefix),
        options=str(CONFIG.third_party_tools.star.align.command_line_arguments),
        prefix=OUTPUT_DIR + "/aligned/star/{sample}_",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/star/{sample}_Aligned.sortedByCoord.out.bam"),
        bam2=temp(
            OUTPUT_DIR + "/aligned/star/{sample}_Aligned.toTranscriptome.out.bam"
        ),
        log_out=temp(OUTPUT_DIR + "/aligned/star/{sample}_Log.final.out"),
    threads: CONFIG.third_party_tools.star.align.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=35, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/align/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}.tsv",
    message: "Aligning reads for sample {wildcards.sample} using STAR",
    shell: """
    STAR \
    --genomeDir {params.index} \
    --readFilesIn {input.fq1} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN {threads} \
    --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} \
    --outFileNamePrefix {params.prefix} \
    {params.options} \
    > {log} 2>&1
    """


rule rename_aligned:
    input:
        bam=OUTPUT_DIR + "/aligned/star/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/raw/{sample}.bam"),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/rename_aligned/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/rename_aligned/{sample}.tsv",
    message: "Renaming aligned BAM for sample {wildcards.sample} to standard format",
    shell: "mv {input.bam} {output.bam}"


localrules:
    rename_aligned,
