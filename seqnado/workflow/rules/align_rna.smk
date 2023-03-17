import seqnado.utils as utils


rule align_paired:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
    params:
        index=config["genome"]["indicies"],
        options=utils.check_options(config["star"]["options"]),
    output:
        bam=temp("seqnado_output/aligned/star/{sample}Aligned.sortedByCoord.out.bam"),
    threads: config["star"]["threads"]
    resources:
        mem_mb=(32000 // config["star"]["threads"]),
    log:
        "seqnado_output/logs/align/{sample}.log",
    shell:
        """
        STAR \
        --genomeDir \
        {params.index} \
        --readFilesIn \
        {input.fq1} \
        {input.fq2} \
        --readFilesCommand \
        zcat \
        --outSAMtype \
        BAM SortedByCoordinate \
        --runThreadN \
        {threads} \
        --outFileNamePrefix \
        aligned/star/{wildcards.sample} \
        {params.options} \
        > {log} 2>&1
        """


rule rename_aligned:
    input:
        bam=rules.align_paired.output.bam,
    output:
        bam="seqnado_output/aligned/sorted/{sample}.bam",
    shell:
        "mv {input.bam} {output.bam}"

localrules: rename_aligned