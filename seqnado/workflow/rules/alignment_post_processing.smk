import ngs_pipeline.utils as utils


rule sort_bam:
    input:
        bam="aligned/{sample}.bam",
    output:
        sentinel=touch("flags/{sample}.sorted.sentinel"),
    threads: 8
    shell:
        """samtools sort {input.bam} -@ {threads} -o {input.bam}_sorted &&
           mv {input.bam}_sorted {input.bam}
        """


rule index_bam:
    input:
        bam="aligned/{sample}.bam",
    output:
        index="aligned/{sample}.bam.bai",
    threads: 1
    shell:
        "samtools index -@ {threads} -b {input.bam}"


rule remove_blacklisted_regions:
    input:
        bam="aligned_and_filtered/{sample}.bam",
    output:
        sentinel=touch("flags/{sample}.blacklist.sentinel"),
    threads: 1
    params:
        blacklist=config["blacklist"],
    log:
        "logs/blacklist/{sample}.log",
    script:
        "../scripts/remove_blacklist.py"


rule shift_atac_alignments:
    input:
        bam="aligned_and_filtered/{sample}.bam",
    output:
        sentinel=touch("flags/{sample}.shifted.sentinel"),
    params:
        options=None,
    threads: 1
    log:
        "logs/duplicate_removal/deeptools/{sample}.log",
    run:
        if config.get("shift_atac_reads"):

            cmd = f"""
                          rsbamtk shift -b {input.bam} -o {input.bam}.tmp &&
                          samtools sort {input.bam}.tmp -@ {threads} -o {input.bam} &&
                          samtools index {input.bam}
                          """

        else:
            cmd = f"""echo "Will not shift reads" > {log}"""
            shell(cmd)


rule mark_filtering_complete:
    input:
        sentinel="flags/{sample}.blacklist.sentinel",
        sentinel2="flags/{sample}.shifted.sentinel",
    output:
        sentinel=touch("flags/{sample}.filtering.complete.sentinel"),
    log:
        "logs/filtering/{sample}.log",
    shell:
        """echo "Filtering complete" > {log}"""
