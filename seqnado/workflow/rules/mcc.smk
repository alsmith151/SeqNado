

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
        samtools sort -@ 4 -o {output.bam}.sorted {output.bam} 2>> {log} &&
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
        bam=temp("seqnado_output/mcc/replicates/{sample}/{sample}_unsorted.bam"),
    params:
        output_dir="seqnado_output/mcc/{sample}/reporters/raw/",
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/identify_viewpoint_reads/{sample}.log",
    container: None
    shell:
        """
        mccnado annotate-bam-file {input.bam} {output.bam} > {log} 2>&1
        """

use rule sort_bam as sort_bam_viewpoints with:
    input:
        bam="seqnado_output/mcc/replicates/{sample}/{sample}_unsorted.bam",
    output:
        bam="seqnado_output/mcc/replicates/{sample}/{sample}.bam",
    log:
        "seqnado_output/logs/sort_bam_viewpoints/{sample}.log",

use rule index_bam as index_bam_viewpoints with:
    input:
        bam="seqnado_output/mcc/replicates/{sample}/{sample}.bam",
    output:
        bai="seqnado_output/mcc/replicates/{sample}/{sample}.bam.bai",


rule extract_ligation_stats:
    input:
        bam="seqnado_output/mcc/replicates/{sample}/{sample}.bam",
    output:
        stats="seqnado_output/resources/{sample}_ligation_stats.json"
    container: None
    shell:
        """
        mccnado extract-ligation-stats {input.bam} {output.stats} 
        """


def get_n_cis_scaling_factor(wc):
    import json 
    from pathlib import Path

    # Can either extract the stats for the sample or the group
    if hasattr(wc, "group"):
        stats_file = f"seqnado_output/resources/{wc.group}_ligation_stats.json"
    else:
        # If not, then use the sample
        stats_file = f"seqnado_output/resources/{wc.sample}_ligation_stats.json"
    
    # Create Path object and ensure the file exists
    stats_path = Path(stats_file)
    if not stats_path.exists():
        raise FileNotFoundError(f"Stats file not found: {stats_file}")
        
    with open(stats_path, 'r') as r:
        stats = json.load(r)
    
    # Check if viewpoint_group exists in stats
    if wc.viewpoint_group not in stats:
        raise KeyError(f"Viewpoint group '{wc.viewpoint_group}' not found in stats file")
    
    # Check if required keys exist in the viewpoint group stats
    required_keys = ['n_cis', 'n_total']
    missing_keys = [key for key in required_keys if key not in stats[wc.viewpoint_group]]
    if missing_keys:
        raise KeyError(f"Missing required keys in stats: {', '.join(missing_keys)}")
    
    # Avoid division by zero
    if stats[wc.viewpoint_group]['n_total'] == 0:
        return 0
        
    return (stats[wc.viewpoint_group]['n_cis'] / stats[wc.viewpoint_group]['n_total']) * 1e6


rule make_bigwigs_mcc_replicates:
    input:
        bam="seqnado_output/mcc/replicates/{sample}/{sample}.bam",
        bai="seqnado_output/mcc/replicates/{sample}/{sample}.bam.bai",
        excluded_regions="seqnado_output/resources/exclusion_regions.bed",
        cis_or_trans_stats="seqnado_output/resources/{sample}_ligation_stats.json",
    output:
        bigwig="seqnado_output/mcc/replicates/{sample}/bigwigs/{viewpoint_group}.bigWig"
    log:
        "seqnado_output/logs/bigwig/{sample}_{viewpoint_group}.log",
    params:
        bin_size=10,
        scale_factor=lambda wc: get_n_cis_scaling_factor(wc),
    container: None
    shell:
        """
        bamnado \
        bam-coverage \
        -b {input.bam} \
        -o {output.bigwig} \
        --bin-size {params.bin_size} \
        --scale-factor {params.scale_factor} \
        --blacklisted-locations {input.excluded_regions} \
        --min-mapq 0 \
        --read-group {wildcards.viewpoint_group} > {log} 2>&1
        """
        
def get_mcc_bam_files_for_merge(wildcards):
    from seqnado.design import NormGroups
    norm_groups = NormGroups.from_design(DESIGN, subset_column="merge")

    sample_names = norm_groups.get_grouped_samples(wildcards.group)
    bam_files = [
        f"seqnado_output/mcc/replicates/{sample}/{sample}.bam" for sample in sample_names
    ]
    return bam_files


rule merge_mcc_bams:
    input:
        bams=get_mcc_bam_files_for_merge,
    output:
        "seqnado_output/mcc/{group}/{group}.bam",
    threads: config["samtools"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/merge_bam/{group}.log",
    shell:
        """
        samtools merge {output} {input} -@ {threads}
        """


use rule index_bam as index_bam_merged with:
    input:
        bam="seqnado_output/mcc/{group}/{group}.bam",
    output:
        bai="seqnado_output/mcc/{group}/{group}.bam.bai",
    log:
        "seqnado_output/logs/index_bam_merged/{group}.log",


use rule extract_ligation_stats as extract_ligation_stats_merged with:
    input:
        bam="seqnado_output/mcc/{group}/{group}.bam",
    output:
        stats="seqnado_output/resources/{group}_ligation_stats.json",
    log:
        "seqnado_output/logs/extract_ligation_stats_merged/{group}.log",


use rule make_bigwigs_mcc_replicates as make_bigwigs_mcc_grouped with:
    input:
        bam="seqnado_output/mcc/{group}/{group}.bam",
        bai="seqnado_output/mcc/{group}/{group}.bam.bai",
        excluded_regions="seqnado_output/resources/exclusion_regions.bed",
        cis_or_trans_stats="seqnado_output/resources/{group}_ligation_stats.json",
    output:
        bigwig="seqnado_output/mcc/{group}/bigwigs/{viewpoint_group}.bigWig"
    params:
        bin_size=10,
        scale_factor=lambda wc: get_n_cis_scaling_factor(wc),
    log:
        "seqnado_output/logs/bigwig/{group}_{viewpoint_group}.log",
    container: None
        


rule identify_ligation_junctions:
    input:
        bam="seqnado_output/mcc/{group}/{group}.bam",
        bai="seqnado_output/mcc/{group}/{group}.bam.bai",
    output:
        pairs=temp(expand("seqnado_output/mcc/{{group}}/ligation_junctions/raw/{viewpoint}.pairs", viewpoint=GROUPED_VIEWPOINT_OLIGOS)),
    log:
        "seqnado_output/logs/ligation_junctions/{group}.log",
    threads: 1
    resources:
        mem="1GB",
    container: None,
    params:
        outdir="seqnado_output/mcc/{group}/ligation_junctions/raw/",
    shell:
        """
        mccnado identify-ligation-junctions \
        {input.bam} \
        {params.outdir}
        """


rule sort_ligation_junctions:
    input:
        pairs="seqnado_output/mcc/{group}/ligation_junctions/raw/{viewpoint}.pairs",
    output:
        pairs=temp("seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.pairs"),
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/sort_ligation_junctions/{group}_{viewpoint}.log",
    shell:
        """
        sort -k2,2 -k4,4 -k3,3n -k5,5n {input.pairs} > {output.pairs}
        """

rule bgzip_pairs:
    input:
        pairs="seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.pairs",
    output:
        pairs="seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.pairs.gz",
    log:
        "seqnado_output/logs/bgzip_pairs/{group}_{viewpoint}.log",
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
        pairs="seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.pairs.gz",
        bins="seqnado_output/resources/genomic_bins.bed",
    output:
        cooler=temp("seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.cool"),
    log:
        "seqnado_output/logs/make_cooler/{group}_{viewpoint}.log",
    params:
        resolution=config.get("resolution", 100),
        genome=config["genome"]["name"],
    container: "library://asmith151/seqnado/seqnado_mcc:latest"
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
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
        cooler="seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.cool",
    output:
        cooler=temp("seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.mcool"),
    log:
        "seqnado_output/logs/zoomify_cooler/{group}_{viewpoint}.log",
    params:
        resolutions=",".join([str(r) for r in config.get("resolutions", [100, 1000, 10000])]),
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
    container: "library://asmith151/seqnado/seqnado_mcc:latest"
    shell:
        """
        cooler zoomify {input.cooler} -r {params.resolutions} -o {output.cooler} > {log} 2>&1
        """


rule aggregate_coolers:
    input:
        mcools=expand("seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.mcool", group=SAMPLE_GROUPS, viewpoint=GROUPED_VIEWPOINT_OLIGOS),
    output:
        mcool="seqnado_output/mcc/{group}/{group}.mcool",
    log:
        "seqnado_output/logs/{group}_aggregate_coolers.log",
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
    container: None
    shell:
        """
        mccnado combine-ligation-junction-coolers \
        {input.mcools} \
        {output.mcool}
        """


rule call_mcc_peaks: # TODO: ensure that we're using the GPU queue
    input:
        bigwig="seqnado_output/mcc/{group}/bigwigs/{viewpoint_group}.bigWig",
    output:
        peaks="seqnado_output/mcc/{group}/peaks/{viewpoint_group}.bed",
    log:
        "seqnado_output/logs/call_mcc_peaks/{group}_{viewpoint_group}.log",
    params:
        options=check_options(config["lanceotron_mcc"]["options"]),
    container: None
    threads: 2
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        gpu=1,
    shell:
        """
        lanceotron-mcc \
        call-mcc-peaks \
        --bigwig {input.bigwig} \
        --outfile {output.peaks} \
        --n-jobs {threads} \
        {params.options} > {log} 2>&1
        """





ruleorder:
    combine_genome_mapped_reads > align_paired




