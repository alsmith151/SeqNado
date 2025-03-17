

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
        bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output.fasta} -nameOnly 2> {log} &&
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
    container: None
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
    container: None
    run:
        from mcc import mcc        
        
        mcc.split_viewpoint_reads(input.bam, output.fq)


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

checkpoint identify_viewpoint_reads:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        bams=directory("seqnado_output/mcc/{sample}/reporters/raw/"),
    params:
        output_dir=lambda wc, output: str(pathlib.Path(output.bams).parent.absolute()),
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/split_genomic_reads/{sample}.log",
    container: None
    run:
        from mcc import mcc
        import pathlib

        outdir = pathlib.Path(params.output_dir)
        outdir.mkdir(exist_ok=True, parents=True)
        mcc.split_genomic_reads(input.bam, params.output_dir)


use rule sort_bam as sort_bam_viewpoints with:
    input:
        bam="seqnado_output/mcc/{sample}/reporters/raw/{sample}.{viewpoint}.bam",
    output:
        bam=temp("seqnado_output/mcc/{sample}/reporters/sorted/{sample}.{viewpoint}.bam"),
    log:
        "seqnado_output/logs/sort_bam/{sample}_{viewpoint}.log",

use rule index_bam as index_bam_viewpoints with:
    input:
        bam="seqnado_output/mcc/{sample}/reporters/sorted/{sample}.{viewpoint}.bam",
    output:
        bai=temp("seqnado_output/mcc/{sample}/reporters/sorted/{sample}.{viewpoint}.bam.bai"),

use rule move_bam_to_final_location as move_bam_viewpoints with:
    input:
        bam="seqnado_output/mcc/{sample}/reporters/sorted/{sample}.{viewpoint}.bam",
    output:
        bam="seqnado_output/mcc/{sample}/reporters/{sample}.{viewpoint}.bam",
        bai="seqnado_output/mcc/{sample}/reporters/{sample}.{viewpoint}.bam.bai",
    log:
        "seqnado_output/logs/move_bam/{sample}_{viewpoint}.log",

use rule deeptools_make_bigwigs as deeptools_make_bigwigs_mcc_replicates with:
    input:
        bam="seqnado_output/mcc/{sample}/reporters/{sample}.{viewpoint}.bam",
        bai="seqnado_output/mcc/{sample}/reporters/{sample}.{viewpoint}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/unscaled/{sample}/{viewpoint}.bigWig",
    log:
        "seqnado_output/logs/deeptools_bigwig/{sample}_{viewpoint}.log",


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
    container: None
    run:
        from mcc import mcc
        import pathlib

        outdir = pathlib.Path(output.pairs)
        outdir.mkdir(exist_ok=True, parents=True)
        mcc.identify_ligation_junctions(input.bam, output.pairs)


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
        None
    shell:
        """
        cooler makebins {input.chrom_sizes} {params.bin_size} > {output.bed}
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
    container: None
    shell:
        """
        cooler cload pairs \
        {input.bins} \
        {input.pairs} \
        {output.cooler} \
        --assembly {params.genome} \
        -c1 2 -p1 3 -c2 4 -p2 5  > 2>$1 {log}
        """

rule zoomify_cooler:
    input:
        cooler="seqnado_output/mcc/{sample}/raw/{viewpoint}.cool",
    output:
        cooler="seqnado_output/mcc/{sample}/{viewpoint}.mcool",
    log:
        "seqnado_output/logs/zoomify_cooler/{sample}_{viewpoint}.log",
    params:
        resolutions=config.get("resolutions", [100, 1000, 10000]),
    container: None
    shell:
        """
        cooler zoomify {input.cooler} {output.cooler} {params.resolutions} > {log}
        """





ruleorder:
    combine_genome_mapped_reads > align_paired




