import re


rule deeptools_make_bigwigs_rna_plus:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/{sample}_plus.bigWig",
    threads: config["deeptools"]["threads"]
    resources:
        mem_mb=500
        time='02:00:00',
    log:
        "seqnado_output/logs/pileups/deeptools/{sample}_plus.log",
    shell:
        """
        bamCoverage -p {threads} --filterRNAstrand forward -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """


rule deeptools_make_bigwigs_rna_minus:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/{sample}_minus.bigWig",
    threads: config["deeptools"]["threads"]
    resources:
        mem_mb=500
        time='02:00:00',
    log:
        "seqnado_output/logs/pileups/deeptools/{sample}_minus.log",
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --filterRNAstrand reverse -p {threads} --scaleFactor -1 > {log} 2>&1
        """
