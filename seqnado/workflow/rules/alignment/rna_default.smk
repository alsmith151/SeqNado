from seqnado.helpers import define_memory_requested, define_time_requested

rule align_paired:
    input:
        fq1=OUTPUT_DIR + "/trimmed/{sample}_1.fastq.gz",
        fq2=OUTPUT_DIR + "/trimmed/{sample}_2.fastq.gz",
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.star.command_line_arguments),
        prefix=OUTPUT_DIR + "/aligned/star/{sample}_",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/star/{sample}_Aligned.sortedByCoord.out.bam"),
        bam2=temp(
            OUTPUT_DIR + "/aligned/star/{sample}_Aligned.toTranscriptome.out.bam"
        ),
        log_out=temp(OUTPUT_DIR + "/aligned/star/{sample}_Log.final.out"),
    threads: config["star"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=35, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/align/{sample}.log",
    shell:
        """
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


rule rename_aligned:
    input:
        bam=rules.align_paired.output.bam,
    output:
        bam=temp(OUTPUT_DIR + "/aligned/raw/{sample}.bam"),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        "mv {input.bam} {output.bam}"


localrules:
    rename_aligned,
