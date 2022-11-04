import ngs_pipeline.utils as utils



rule align_paired:
    input:
        fq1="trimmed/{sample}_1.fastq.gz",
        fq2="trimmed/{sample}_2.fastq.gz",
    params:
        index=config["genome"]["indicies"],
        options=utils.check_options(config["star"]["options"]),
    output:
        bam="aligned/{sample}Aligned.sortedByCoord.out.bam",
    threads: config["star"]["threads"]
    log:
        "logs/align/{sample}.log",
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
        aligned/{wildcards.sample} \
        {params.options} 
        """


rule rename_aligned:
    input:
        bam = "aligned/{sample}Aligned.sortedByCoord.out.bam",
    output:
        bam = "aligned/{sample}.bam"
    shell:
        "mv {input.bam} {output.bam}"