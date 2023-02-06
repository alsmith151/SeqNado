import re


rule deeptools_make_bigwigs_rna_plus:
    input:
        bam="aligned_and_filtered/{sample}.bam",
        bai="aligned_and_filtered/{sample}.bam.bai",
        filtering="flags/{sample}.filtering.complete.sentinel",
    output:
        bigwig="bigwigs/deeptools/{sample}_plus.bigWig",
    threads: config["deeptools"]["threads"]
    log:
        "logs/pileups/deeptools/{sample}_plus.log",
    shell:
        """
        bamCoverage -p {threads} --filterRNAstrand forward -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """


rule deeptools_make_bigwigs_rna_minus:
    input:
        bam="aligned_and_filtered/{sample}.bam",
        bai="aligned_and_filtered/{sample}.bam.bai",
        filtering="flags/{sample}.filtering.complete.sentinel",
    output:
        bigwig="bigwigs/deeptools/{sample}_minus.bigWig",
    threads: config["deeptools"]["threads"]
    log:
        "logs/pileups/deeptools/{sample}_minus.log",
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --filterRNAstrand reverse -p {threads} > {log} 2>&1
        """
