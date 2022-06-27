rule mark_duplicates:
    input:
        bam="aligned/{sample}.bam",
    output:
        bam="aligned_and_filtered/{sample}.bam",
        metrics="aligned_and_filtered/{sample}.metrics.txt",
        log = "logs/duplicate_removal/picard/{sample}.log",
    log:
        "logs/duplicate_removal/picard/{sample}.log",
    params:
        extra = "--REMOVE_DUPLICATES true --CREATE_INDEX true",
    resources:
        mem_mb=1024 * 4,
    threads:
        4,
    shell:
        """
        picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} {params.extra} > {log} 2>&1
        """


