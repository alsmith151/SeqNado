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
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        """
        mccnado annotate-bam-file {input.bam} {output.bam} > {log} 2>&1
        """

use rule sort_bam as sort_bam_viewpoints with:
    input:
        bam="seqnado_output/mcc/replicates/{sample}/{sample}_unsorted.bam",
    output:
        bam="seqnado_output/mcc/replicates/{sample}/{sample}.bam",
        read_log="seqnado_output/mcc/replicates/{sample}/{sample}_read_log.txt",
    log:
        "seqnado_output/logs/sort_bam_viewpoints/{sample}.log",

use rule index_bam as index_bam_viewpoints with:
    input:
        bam="seqnado_output/mcc/replicates/{sample}/{sample}.bam",
    output:
        bai="seqnado_output/mcc/replicates/{sample}/{sample}.bam.bai",

# ruleorder:
#     combine_genome_mapped_reads > align_paired




