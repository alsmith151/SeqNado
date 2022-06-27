
rule deeptools_filter_alignments:
    input:
        bam = "aligned/{sample}.bam",
        index = "aligned/{sample}.bam.bai"
    output:
        bam = "aligned_and_filtered/{sample}.bam",
        index = "aligned_and_filtered/{sample}.bam.bai",
        log = "logs/duplicate_removal/deeptools/{sample}.log",
    params:
        deduplicate = "--ignoreDuplicates",
        options = config["deeptools"]["alignmentsieve"],
    threads:
        8
    log:
        "logs/duplicate_removal/deeptools/{sample}.log",
    shell:
        """
        alignmentSieve -b {input.bam} -o {output.bam} -p {threads} {params.deduplicate} {params.options} > {log} 2>&1 &&
        samtools sort -o {output.bam}.tmp {output.bam} -@ {threads} &&
        mv {output.bam}.tmp {output.bam} &&
        samtools index {output.bam} &&
        rm -f {output.bam}.tmp
        """