from seqnado.helpers import check_options

rule align_paired:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
    params:
        index=config["genome"]["indices"],
        options=check_options(config["star"]["options"]),
        prefix="seqnado_output/aligned/star/{sample}_",
    output:
        bam=temp("seqnado_output/aligned/star/{sample}_Aligned.sortedByCoord.out.bam"),
        bam2=temp(
            "seqnado_output/aligned/star/{sample}_Aligned.toTranscriptome.out.bam"
        ),
    threads: config["star"]["threads"]
    resources:
        mem="35GB",
        runtime="6h",
    log:
        "seqnado_output/logs/align/{sample}.log",
    shell:
        """
        STAR \
        --genomeDir {params.index} \
        --readFilesIn {input.fq1} {input.fq2} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --runThreadN {threads} \
        --outFileNamePrefix {params.prefix} {params.options} > {log} 2>&1
        """


rule rename_aligned:
    input:
        bam=rules.align_paired.output.bam,
    output:
        bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
    shell:
        "mv {input.bam} {output.bam}"


localrules:
    rename_aligned,
