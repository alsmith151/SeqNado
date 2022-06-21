rule mark_duplicates:
    input:
        bam="aligned/{sample}.bam",
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam="aligned_and_filtered/{sample}.bam",
        metrics="aligned_and_filtered/{sample}.metrics.txt",
        log = "logs/duplicate_removal/picard/{sample}.log",
    log:
        "logs/duplicate_removal/picard/{sample}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024 * 4,
    threads:
        4,
    shell:
        """
        picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} {params.extra} > {log} 2>&1
        """
        
    
rule index_bam:
    input:
        bam="aligned_and_filtered/{sample}.bam",
    output:
        index="aligned_and_filtered/{sample}.bam.bai",
    threads:
        1
    shell:
        "samtools index {input.bam} -@ {threads}"

