import re
from seqnado.helpers import  define_time_requested, define_memory_requested
from seqnado.config.third_party_tools import CommandLineArguments

def format_deeptools_options(wildcards, options: CommandLineArguments) -> str:
    """
    Format the command line options for deeptools based on the input files and parameters.

    Mainly this removes the extend reads option if single ended
    """
    search_term = f'{wildcards.sample}'
    is_paired = INPUT_FILES.is_paired_end(search_term)
    if not is_paired:
        options = CommandLineArguments(value=options, exclude={"--extendReads", "-e", "--samFlagInclude 3"})

    return str(options)


rule homer_make_tag_directory:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
    output:
        homer_tag_directory=directory("seqnado_output/tag_dirs/{sample}"),
    params:
        options=check_options(config["homer"]["maketagdirectory"]),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
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
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
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
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    threads:
        config["deeptools"]["threads"]
    log:
        "seqnado_output/logs/pileups/deeptools/unscaled/{sample}.log",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
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
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/pileups/deeptools/unscaled/{sample}_plus.log",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
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
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/pileups/deeptools/unscaled/{sample}_minus.log",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        """
        bamCoverage {params.options} -p {threads} -b {input.bam} -o {output.bigwig} --filterRNAstrand reverse --scaleFactor -1 > {log} 2>&1
        """

rule bamnado_bam_coverage:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/bamnado/unscaled/{sample}.bigWig",
    params:
        options=check_options(config["bamnado"]["bamcoverage"]),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    threads: config["bamnado"]["threads"],
    log:
        "seqnado_output/logs/pileups/bamnado/{sample}.log",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """

rule bamnado_bam_coverage_rna_plus:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/bamnado/{sample}_plus.bigWig",
    params:
        options=check_options(config["bamnado"]["bamcoverage"]),
    threads: config["bamnado"]["threads"],
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/pileups/bamnado/{sample}_plus.log",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand forward > {log} 2>&1
        """

rule bamnado_bam_coverage_rna_minus:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/bamnado/{sample}_minus.bigWig",
    params:
        options=check_options(config["bamnado"]["bamcoverage"]),
    threads: config["bamnado"]["threads"],
    resources: 
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/pileups/bamnado/{sample}_minus.log",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand reverse --scale-factor -1 > {log} 2>&1
        """

rule fragment_bedgraph:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        filtered=temp("seqnado_output/bedgraphs/{sample}.filtered.bam"),
        sort=temp("seqnado_output/bedgraphs/{sample}.sorted.bam"),
        bed=temp("seqnado_output/bedgraphs/{sample}.bed"),
        bed_log=temp("seqnado_output/logs/bedgraphs/{sample}_bamtobed.log"),
        fragments=temp("seqnado_output/bedgraphs/{sample}.fragments.bed"),
        bdg=temp("seqnado_output/bedgraphs/{sample}.bedGraph"),
    params:
        genome=config['genome']['chromosome_sizes'],
    threads: 16
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=12, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),    
    log:
        "seqnado_output/logs/bedgraphs/{sample}.log",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:"""
        samtools view -@ {threads} -q 30 -f 2 -h {input.bam} | grep -v chrM > {output.filtered} 2> {log}
        samtools sort -@ {threads} -m 900M -o {output.sort} -T {output.sort}.tmp {output.filtered} 2>> {log}
        bedtools bamtobed -bedpe -i {output.sort} > {output.bed} 2>> {output.bed_log}
        awk '$1==$4 && $6-$2 < 1000' {output.bed} > {output.fragments}.temp 2>> {log}
        awk 'BEGIN {{OFS="\t"}} {{print $1, $2, $6}}' {output.fragments}.temp | sort -k1,1 -k2,2n -k3,3n > {output.fragments} 2>> {log}
        bedtools genomecov -bg -i {output.fragments} -g {params.genome} > {output.bdg} 2>> {log}
        rm seqnado_output/bedgraphs/{wildcards.sample}.fragments.bed.temp
        """

ruleorder: deeptools_make_bigwigs_rna_plus > deeptools_make_bigwigs_rna_minus > deeptools_make_bigwigs
