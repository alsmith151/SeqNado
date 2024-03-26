import re
from seqnado.helpers import check_options

def format_deeptools_options(wildcards, options):
    is_paired = DESIGN.query(wildcards.sample).is_paired

    if not is_paired:
        options = re.sub(r"--extendReads", "", options)
        options = re.sub(r"-e", "", options)
    if not options:
        return ""
    else:
        return options


rule homer_make_tag_directory:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
    output:
        homer_tag_directory=directory("seqnado_output/tag_dirs/{sample}"),
    params:
        options=check_options(config["homer"]["maketagdirectory"]),
    resources:
        mem="4GB",
        runtime="2h",
    log:
        "seqnado_output/logs/homer/maketagdirectory_{sample}.log",
    shell:
        """makeTagDirectory {output.homer_tag_directory} {input.bam} {params.options} > {log} 2>&1"""


rule homer_make_bigwigs:
    input:
        homer_tag_directory="seqnado_output/tag_dirs/{sample}",
    output:
        homer_bigwig="seqnado_output/bigwigs/homer/unscaled/{sample}.bigWig",
    log:
        "seqnado_output/logs/homer/makebigwigs_{sample}.log",
    resources:
        mem="4GB",
        runtime="2h",
    params:
        genome_name=config["genome"]["name"],
        genome_chrom_sizes=config["genome"]["chromosome_sizes"],
        options=check_options(config["homer"]["makebigwig"]),
        outdir="seqnado_output/bigwigs/homer/",
        temp_bw=lambda wc, output: output.homer_bigwig.replace(
            ".bigWig", ".ucsc.bigWig"
        ),
    shell:
        """makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir {params.outdir} {params.options} > {log} 2>&1 &&
           mv {params.outdir}/{wildcards.sample}.ucsc.bigWig {output.homer_bigwig}
        """


rule deeptools_make_bigwigs:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/unscaled/{sample}.bigWig",
    params:
        options=lambda wildcards: format_deeptools_options(wildcards, config["deeptools"]["bamcoverage"]),
    resources:
        mem="2GB",
        runtime="4h",
    threads: 
        config["deeptools"]["threads"]
    log:
        "seqnado_output/logs/pileups/deeptools/{sample}.log",
    shell:
        """
        bamCoverage {params.options} -p {threads} -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """


rule deeptools_make_bigwigs_rna_plus:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/unscaled/{sample}_plus.bigWig",
    params:
        options=check_options(config["deeptools"]["bamcoverage"]),
    threads: config["deeptools"]["threads"]
    resources:
        mem="2GB",
        runtime="4h",
    log:
        "seqnado_output/logs/pileups/deeptools/{sample}_plus.log",
    shell:
        """
        bamCoverage {params.options} -p {threads} --filterRNAstrand forward -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """


rule deeptools_make_bigwigs_rna_minus:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/unscaled/{sample}_minus.bigWig",
    params:
        options=check_options(config["deeptools"]["bamcoverage"]),
    threads: config["deeptools"]["threads"]
    resources:
        mem="2GB",
        runtime="4h",
    log:
        "seqnado_output/logs/pileups/deeptools/{sample}_minus.log",
    shell:
        """
        bamCoverage {params.options} -p {threads} -b {input.bam} -o {output.bigwig} --filterRNAstrand reverse --scaleFactor -1 > {log} 2>&1
        """


rule fragment_bedgraph:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        bedgraph="seqnado_output/bedgraphs/{sample}.bedGraph",
    params:
        genome=config['genome']['chromosome_sizes']
        outdir="seqnado_output/bedgraphs/",
    shell:
        """
        bedtools bamtobed -bedpe -i {input.bam} > {params.outdir}/{wildcards.sample}.bed
        awk '$1==$4 && $6-$2 < 1000 {print $0}' {params.outdir}/{wildcards.sample}.bed > {params.outdir}/{wildcards.sample}_clean.bed
        cut -f 1,2,6 {params.outdir}/{wildcards.sample}_clean.bed | sort -k1,1 -k2,2n -k3,3n > {params.outdir}/{wildcards.sample}_fragments.bed
        bedtools genomecov -bg -i {params.outdir}/{wildcards.sample}_fragments.bed -g {params.genome} > {output.bedgraph}
        
        rm {params.outdir}/{wildcards.sample}.bed
        rm {params.outdir}/{wildcards.sample}_clean.bed
        rm {params.outdir}/{wildcards.sample}_fragments.bed
        """"

ruleorder: deeptools_make_bigwigs_rna_plus > deeptools_make_bigwigs_rna_minus > deeptools_make_bigwigs
