import re
import seqnado.utils as utils


rule homer_make_tag_directory:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
    output:
        homer_tag_directory=directory("seqnado_output/tag_dirs/{sample}"),
    params:
        options=utils.check_options(config["homer"]["maketagdirectory"]),
    resources:
        mem_mb=4000
    log:
        "seqnado_output/logs/homer/maketagdirectory_{sample}.log",
    shell:
        """makeTagDirectory {output.homer_tag_directory} {input.bam} {params.options} > {log} 2>&1"""


rule homer_make_bigwigs:
    input:
        homer_tag_directory="seqnado_output/tag_dirs/{sample}",
    output:
        homer_bigwig="seqnado_output/bigwigs/homer/{sample}.bigWig",
    log:
        "seqnado_output/logs/homer/makebigwigs_{sample}.log",
    resources:
        mem_mb=4000
    params:
        genome_name=config["genome"]["name"],
        genome_chrom_sizes=config["genome"]["chromosome_sizes"],
        options=utils.check_options(config["homer"]["makebigwig"]),
        outdir="seqnado_output/bigwigs/homer/",
        temp_bw=lambda wc, output: output.homer_bigwig.replace(".bigWig", ".ucsc.bigWig"),
    shell:
        """makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir {params.outdir} {params.options} > {log} 2>&1 &&
           mv {params.outdir}/{wildcards.sample}.ucsc.bigWig {output.homer_bigwig}
        """


rule deeptools_make_bigwigs:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/{sample}.bigWig",
    params:
        options=utils.check_options(config["deeptools"]["bamcoverage"]),
    resources:
        mem_mb=1000
    threads: config["deeptools"]["threads"]
    log:
        "seqnado_output/logs/pileups/deeptools/{sample}.log",
    shell:
        """
        bamCoverage {params.options} -b {input.bam} -o {output.bigwig} -p {threads} > {log} 2>&1
        """
