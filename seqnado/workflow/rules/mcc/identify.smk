def identify_extracted_bam_files(wildcards):
    from pathlib import Path

    checkpoint_output = checkpoints.identify_viewpoint_reads.get(**wildcards)
    outdir = Path(checkpoint_output.output.bams)
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
        outdir = Path(checkpoint_output.output.bams)
        viewpoints = glob_wildcards(str(outdir / "{viewpoint}.bam")).viewpoint
        
        if ii == 0:
            viewpoint_set = set(viewpoints)
        else:
            viewpoint_set = viewpoint_set.intersection(viewpoints)
    return list(viewpoint_set)


rule identify_viewpoint_reads:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
    output:
        bam=temp(OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}_unsorted.bam"),
    params:
        output_dir=OUTPUT_DIR + "/mcc/{sample}/reporters/raw/",
    threads: 1
    resources:
        mem="1GB",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/identify_viewpoint_reads/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/identify_viewpoint_reads/{sample}.tsv",
    message: "Identifying viewpoint reads for sample {wildcards.sample}",
    shell:
        """
        mccnado annotate-bam-file {input.bam} {output.bam} > {log} 2>&1
        """

use rule sort_bam as sort_bam_viewpoints with:
    input:
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}_unsorted.bam",
    output:
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam",
        read_log=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}_read_log.txt",
    log: OUTPUT_DIR + "/logs/sort_bam_viewpoints/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/sort_bam_viewpoints/{sample}.tsv",
    message: "Sorting BAM file for viewpoints for sample {wildcards.sample}",

use rule index_bam as index_bam_viewpoints with:
    input:
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam",
    output:
        bai=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam.bai",
    log: OUTPUT_DIR + "/logs/index_bam_viewpoints/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/index_bam_viewpoints/{sample}.tsv",
    message: "Indexing BAM file for viewpoints for sample {wildcards.sample}",

# ruleorder:
#     combine_genome_mapped_reads > align_paired




