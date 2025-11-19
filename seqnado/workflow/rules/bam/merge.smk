# Combine bam files for merge

def get_bam_files_for_consensus(wildcards):
    """Get BAM files for merging based on sample names."""
    groups = SAMPLE_GROUPINGS.groupings.get(wildcards.group)
    sample_names = groups.get_samples()
    bam_files = [
        OUTPUT_DIR + "/aligned/{sample}.bam" for sample in sample_names
    ]
    return bam_files


rule merge_bams:
    input:
        bams=get_bam_files_for_consensus,
    output:
        temp(OUTPUT_DIR + "/aligned/merged/{group}.bam"),
    threads: CONFIG.third_party_tools.samtools.merge.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/merge_bam/{group}.log",
    shell:"""
    samtools merge {output} {input} -@ {threads}
    """


use rule index_bam as index_consensus_bam with:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
    output:
        bai=temp(OUTPUT_DIR + "/aligned/merged/{group}.bam.bai"),
    threads: 8