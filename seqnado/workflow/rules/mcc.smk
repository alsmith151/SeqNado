

def identify_extracted_bam_files(wildcards):
    import pathlib

    checkpoint_output = checkpoints.identify_viewpoint_reads.get(**wildcards)
    outdir = pathlib.Path(checkpoint_output.output.bams)
    viewpoints = glob_wildcards(str(outdir / "{viewpoint}.bam")).viewpoint
    return expand(str(outdir / "{viewpoint}.bam"), viewpoint=viewpoints)

def redefine_viewpoints(samples):
    """
    Redefine the set of viewpoints to be the intersection of viewpoints across all samples.

    The issue is that some viewpoints may not be present in all samples or may not have enough reads to be considered.

    Parameters
    ----------
    samples : list
        List of samples.
    """

    viewpoint_set = set()
    
    for ii, sample in enumerate(samples):
        checkpoint_output = checkpoints.identify_viewpoint_reads.get(sample=sample)
        outdir = pathlib.Path(checkpoint_output.output.bams)
        viewpoints = glob_wildcards(str(outdir / "{viewpoint}.bam")).viewpoint
        
        if ii == 0:
            viewpoint_set = set(viewpoints)
        else:
            viewpoint_set = viewpoint_set.intersection(viewpoints)
    return list(viewpoint_set)



rule viewpoints_to_fasta:
    input:
        bed=config["viewpoints"],
        genome=config["fasta"],
    output:
        fasta="seqnado_output/resources/viewpoints.fa",
    log:
        "seqnado_output/logs/bed_to_fasta/viewpoints.log",
    
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output.fasta} -name 2> {log} &&
        cat {output.fasta} | sed -E 's/:+/-/g' > {output.fasta}.tmp &&
        mv {output.fasta}.tmp {output.fasta}
        """


rule fasta_index:
    input:
        fasta="seqnado_output/resources/viewpoints.fa",
    output:
        index="seqnado_output/resources/viewpoints.fa.fai",
    log:
        "seqnado_output/logs/bed_to_fasta/index.log",
    shell:
        """
        samtools faidx {input.fasta} -o {output.index}
        """

rule exclusion_regions:
    input:
        bed=config['viewpoints'],
    output:
        bed="seqnado_output/resources/exclusion_regions.bed"
    log:
        "seqnado_output/logs/exclusion_regions.log"
    params:
        genome=config['genome']['chromosome_sizes'],
        exclusion_zone=config.get("exclusion_zone", 500)
    shell:
        """
        bedtools slop -i {input.bed} -g {params.genome}  -b {params.exclusion_zone} > {output.bed}
        """


rule minimap2_to_viewpoints:
    input:
        fq="seqnado_output/flashed/{sample}/{sample}.extendedFrags.fastq.gz",
        viewpoints="seqnado_output/resources/viewpoints.fa",
    output:
        bam=temp("seqnado_output/aligned/aligned_to_viewpoints/{sample}.bam"),
    threads: 4
    resources:
        mem="4GB",
    container: "library://asmith151/seqnado/seqnado_mcc:latest"
    log:
        "seqnado_output/logs/aligned/{sample}.log",
    shell:
        """
        minimap2 -x sr -a -k 8 -w 1 --cs=long {input.viewpoints} {input.fq} 2> {log} |
        samtools view -h -b -o {output.bam} 2>> {log} &&
        samtools sort -o {output.bam}.sorted {output.bam} 2>> {log} &&
        mv {output.bam}.sorted {output.bam} &&
        samtools index {output.bam}
        """ 


rule split_reads_aligned_to_viewpoints:
    input:
        bam="seqnado_output/aligned/aligned_to_viewpoints/{sample}.bam",
    output:
        fq="seqnado_output/mcc/{sample}/{sample}.sliced.fastq.gz",
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/split_reads/{sample}.log",
    container: "library://asmith151/seqnado/seqnado_mcc:latest"
    script:
        "../scripts/mcc_split_reads_aligned_to_viewpoints.py"


use rule align_single as align_mcc_reads_to_genome with:
    input:
        fq1="seqnado_output/mcc/{sample}/{sample}.sliced.fastq.gz", 
    output:
        bam=temp("seqnado_output/aligned/initial_alignment/{sample}.bam"),


rule align_unmapped_reads_to_genome:
    input:
        bam="seqnado_output/aligned/initial_alignment/{sample}.bam",
    output:
        bam=temp("seqnado_output/aligned/second_alignment/{sample}.bam"),
        bai=temp("seqnado_output/aligned/second_alignment/{sample}.bam.bai"),
    threads: config["samtools"]["threads"]
    resources:
        mem="4GB",
    log:
        "seqnado_output/logs/realign/{sample}.log",
    params:
        index=config["genome"]["index"],
        options=check_options(config["bowtie2"]["options"]),

    shell:
        """
        samtools view -b -f 4 {input.bam} | bowtie2 -p {threads} -x {params.index} -b - --very-sensitive-local 2>> {log} |
        samtools view -bS - > {output.bam} &&
        samtools sort -@ {threads} -o {output.bam}_sorted {output.bam} &&
        mv {output.bam}_sorted {output.bam} &&
        samtools index {output.bam}
        """

rule combine_genome_mapped_reads:
    input:
        bam1=rules.align_mcc_reads_to_genome.output.bam,
        bam2=rules.align_unmapped_reads_to_genome.output.bam,
    output:
        bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
    threads: config["samtools"]["threads"]
    resources:
        mem="4GB",
    log:
        "seqnado_output/logs/combine/{sample}.log",
    shell:
        """
        samtools merge -@ {threads} {output.bam} {input.bam1} {input.bam2} &&
        samtools view -F 4 -b {output.bam} > {output.bam}.tmp &&
        mv {output.bam}.tmp {output.bam} &&
        samtools sort -n -@ {threads} -o {output.bam}_sorted {output.bam} &&
        mv {output.bam}_sorted {output.bam}
        """

rule identify_viewpoint_reads:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        bams=expand("seqnado_output/mcc/{{sample}}/reporters/raw/{viewpoint}.bam", viewpoint=VIEWPOINT_OLIGOS),
    params:
        output_dir="seqnado_output/mcc/{sample}/reporters/raw/",
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/split_genomic_reads/{sample}.log",
    container: "library://asmith151/seqnado/seqnado_mcc:latest"
    script:
        "../scripts/mcc_identify_viewpoint_reads.py"


use rule sort_bam as sort_bam_viewpoints with:
    input:
        bam="seqnado_output/mcc/{sample}/reporters/raw/{viewpoint}.bam",
    output:
        bam=temp("seqnado_output/mcc/{sample}/reporters/sorted/{viewpoint}.bam"),
    log:
        "seqnado_output/logs/sort_bam/{sample}_{viewpoint}.log",

use rule index_bam as index_bam_viewpoints with:
    input:
        bam="seqnado_output/mcc/{sample}/reporters/sorted/{viewpoint}.bam",
    output:
        bai=temp("seqnado_output/mcc/{sample}/reporters/sorted/{viewpoint}.bam.bai"),

use rule move_bam_to_final_location as move_bam_viewpoints with:
    input:
        bam="seqnado_output/mcc/{sample}/reporters/sorted/{viewpoint}.bam",
        bai="seqnado_output/mcc/{sample}/reporters/sorted/{viewpoint}.bam.bai",
    output:
        bam="seqnado_output/mcc/{sample}/reporters/{viewpoint}.bam",
        bai="seqnado_output/mcc/{sample}/reporters/{viewpoint}.bam.bai",
    log:
        "seqnado_output/logs/move_bam/{sample}_{viewpoint}.log",

use rule deeptools_make_bigwigs as deeptools_make_bigwigs_mcc_replicates with:
    input:
        bam="seqnado_output/mcc/{sample}/reporters/{viewpoint}.bam",
        bai="seqnado_output/mcc/{sample}/reporters/{viewpoint}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/unscaled/{sample}/{viewpoint}.bigWig",
    log:
        "seqnado_output/logs/deeptools_bigwig/{sample}_{viewpoint}.log",


def define_bigwigs(wc):
    """
    Define the bigwig files to be created for each viewpoint group.
    """
    from collections import defaultdict
    viewpoints_required = defaultdict(list)
    
    # Have a mapping that goes from viewpoint to viewpoint group. Want to collate all the viewpoints in a group.
    for viewpoint, viewpoint_group in VIEWPOINT_TO_GROUPED_VIEWPOINT.items():
        viewpoints_required[viewpoint_group].append(viewpoint)

    # Select the viewpoint group to be used
    viewpoint_group = wc.viewpoint_group
    viewpoints = viewpoints_required[viewpoint_group]

    return expand("seqnado_output/bigwigs/deeptools/unscaled/{sample}/{viewpoint}.bigWig", sample=wc.sample, viewpoint=viewpoints)


rule merge_viewpoint_bigwigs:
    input:
        bigwigs=define_bigwigs,
    output:
        bigwig="seqnado_output/bigwigs/deeptools/grouped_viewpoints/{sample}/{viewpoint_group}.bigWig",
    log:
        "seqnado_output/logs/merge_bigwigs/{sample}_{viewpoint_group}.log",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads:
        config["deeptools"]["threads"],
    params:
        options=lambda wildcards: check_options(config["deeptools"]["bamcoverage"]),
    shell:
        """
        bigwigAverage -b {input.bigwigs} -o {output.bigwig} -p {threads} {params.options} 2> {log}
        """



def get_bigwigs_to_merge(wc):
    from seqnado.design import NormGroups
    norm_groups = NormGroups.from_design(DESIGN, subset_column="merge")

    sample_names = norm_groups.get_grouped_samples(wc.group)

    return expand("seqnado_output/bigwigs/deeptools/grouped_viewpoints/{sample}/{viewpoint_group}.bigWig", sample=sample_names, viewpoint_group=wc.viewpoint_group)

use rule merge_viewpoint_bigwigs as merge_samples_bigwigs with:
    input:
        bigwigs=get_bigwigs_to_merge,
    output:
        bigwig="seqnado_output/bigwigs/deeptools/grouped_samples/{group}_{viewpoint_group}.bigWig",
    log:
        "seqnado_output/logs/merge_samples_bigwigs/{group}_{viewpoint_group}.log",


rule identify_ligation_junctions:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        pairs=expand("seqnado_output/mcc/{{sample}}/ligation_junctions/raw/{viewpoint}.pairs", viewpoint=VIEWPOINT_OLIGOS),
    log:
        "seqnado_output/logs/ligation_junctions/{sample}.log",
    threads: 1
    resources:
        mem="1GB",
    container: "library://asmith151/seqnado/seqnado_mcc:latest"
    params:
        outdir="seqnado_output/mcc/{sample}/ligation_junctions/raw/",
    script:
        "../scripts/mcc_identify_ligation_junctions.py"


rule sort_ligation_junctions:
    input:
        pairs="seqnado_output/mcc/{sample}/ligation_junctions/raw/{viewpoint}.pairs",
    output:
        pairs=temp("seqnado_output/mcc/{sample}/ligation_junctions/{viewpoint}.pairs"),
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/sort_ligation_junctions/{sample}_{viewpoint}.log",
    shell:
        """
        sort -k2,2 -k4,4 -k3,3n -k5,5n {input.pairs} > {output.pairs}
        """

rule bgzip_pairs:
    input:
        pairs="seqnado_output/mcc/{sample}/ligation_junctions/{viewpoint}.pairs",
    output:
        pairs="seqnado_output/mcc/{sample}/ligation_junctions/{viewpoint}.pairs.gz",
    log:
        "seqnado_output/logs/bgzip_pairs/{sample}_{viewpoint}.log",
    shell:
        """
        bgzip -c {input.pairs} > {output.pairs}
        """

rule make_genomic_bins:
    input:
        chrom_sizes=config["genome"]["chromosome_sizes"],
    params:
        bin_size=config["resolution"],
    output:
        bed="seqnado_output/resources/genomic_bins.bed",
    log:
        "seqnado_output/logs/genomic_bins.log",
    container:
        "library://asmith151/seqnado/seqnado_mcc:latest"
    shell:
        """
        cooler makebins {input.chrom_sizes} {params.bin_size} -o {output.bed} > {log} 2>&1
        """


rule make_cooler:
    input:
        pairs="seqnado_output/mcc/{sample}/ligation_junctions/{viewpoint}.pairs.gz",
        bins="seqnado_output/resources/genomic_bins.bed",
    output:
        cooler="seqnado_output/mcc/{sample}/raw/{viewpoint}.cool",
    log:
        "seqnado_output/logs/make_cooler/{sample}_{viewpoint}.log",
    params:
        resolution=config.get("resolution", 100),
        genome=config["genome"]["name"],
    container: "library://asmith151/seqnado/seqnado_mcc:latest"
    shell:
        """
        cooler cload pairs \
        {input.bins} \
        {input.pairs} \
        {output.cooler} \
        --assembly {params.genome} \
        -c1 2 -p1 3 -c2 4 -p2 5 > {log} 2>&1
        """

rule zoomify_cooler:
    input:
        cooler="seqnado_output/mcc/{sample}/raw/{viewpoint}.cool",
    output:
        cooler="seqnado_output/mcc/{sample}/{viewpoint}.mcool",
    log:
        "seqnado_output/logs/zoomify_cooler/{sample}_{viewpoint}.log",
    params:
        resolutions=",".join([str(r) for r in config.get("resolutions", [100, 1000, 10000])]),
    container: "library://asmith151/seqnado/seqnado_mcc:latest"
    shell:
        """
        cooler zoomify {input.cooler} -r {params.resolutions} -o {output.cooler} > {log} 2>&1
        """





ruleorder:
    combine_genome_mapped_reads > align_paired




