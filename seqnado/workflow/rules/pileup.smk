import re
import seqnado.utils as utils


rule homer_make_tag_directory:
    input:
        bam="aligned/{sample}.bam",
    output:
        homer_tag_directory=directory("tag_dirs/{sample}"),
    params:
        options=config["homer"]["maketagdirectory"],
    log:
        "logs/homer/maketagdirectory_{sample}.log",
    shell:
        """makeTagDirectory {output.homer_tag_directory} {input.bam} {params.options} > {log} 2>&1"""


rule homer_make_bigwigs:
    input:
        homer_tag_directory="tag_dirs/{sample}",
    output:
        homer_bigwig="bigwigs/homer/{sample}.bigWig",
    log:
        "logs/homer/makebigwigs_{sample}.log",
    params:
        genome_name=config["genome"]["name"],
        genome_chrom_sizes=config["genome"]["chromosome_sizes"],
        options=config["homer"]["makebigwig"],
        temp_bw=lambda wc, output: output.homer_bigwig.replace(".bigWig", ".ucsc.bigWig"),
    shell:
        """makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir bigwigs/homer/ {params.options} > {log} 2>&1 &&
           mv {params.temp_bw} {output.homer_bigwig}
        """


rule deeptools_make_bigwigs:
    input:
        bam="aligned/{sample}.bam",
        bai="aligned/{sample}.bam.bai",
    output:
        bigwig="bigwigs/deeptools/{sample}.bigWig",
    params:
        options=config["deeptools"]["bamcoverage"],
    threads: config["deeptools"]["threads"]
    log:
        "logs/pileups/deeptools/{sample}.log",
    shell:
        """
        bamCoverage {params.options} -b {input.bam} -o {output.bigwig} -p {threads} > {log} 2>&1
        """
