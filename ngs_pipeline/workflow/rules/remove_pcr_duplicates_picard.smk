rule mark_duplicates:
    input:
        bam="aligned/{sample}.bam",
        index="aligned/{sample}.bam.bai",
    output:
        bam="aligned_and_filtered/{sample}.bam",
        metrics="aligned_and_filtered/{sample}.metrics.txt",
    log:
        "logs/duplicate_removal/picard/{sample}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    resources:
        mem_mb=1024 * 4,
    threads: 4
    shell:
        """
        picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} {params.extra} > {log} 2>&1
        """


rule reindex_bam:
    input:
        bam="aligned_and_filtered/{sample}.bam",
        filtering="logs/blacklist/{sample}.log",
    output:
        index="aligned_and_filtered/{sample}.bam.bai",
    threads: 1
    shell:
        "samtools index {input.bam} -@ {threads}"
