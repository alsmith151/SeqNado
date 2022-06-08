rule deeptools_filter_alignments:
    input:
        bam = "aligned_and_filtered/{sample}.bam",
    output:
        log = "logs/shift_alignments/{sample}.log",
    params:
        options = config["deeptools"]["alignmentsieve"].replace("--ignoreDuplicates", ""),
    threads:
        8
    log:
        "logs/duplicate_removal/deeptools/{sample}.log",
    shell:
        """
        alignmentSieve -b {input.bam} -o {output.bam} -p {threads} {params.deduplicate} {params.options} --ATACshift > {log} 2>&1 &&
        samtools sort -o {output.bam}.tmp {output.bam} -@ {threads} &&
        mv {output.bam}.tmp {output.bam} &&
        samtools index {output.bam} &&
        rm -f {output.bam}.tmp
        """