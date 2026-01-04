"""
Rules for comparing MCC data across conditions by merging BAM files.
"""

rule merge_condition_mcc_bams:
    input:
        bams=get_mcc_bam_files_for_merge,
    output:
        OUTPUT_DIR + "/mcc/{condition}/{condition}.bam",
    threads: CONFIG.third_party_tools.samtools.merge.threads
    wildcard_constraints:
        condition="|".join(SAMPLE_GROUPINGS.get_grouping('condition').group_names),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/merge_bam/{condition}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/merge_bam/{condition}.tsv",
    message: "Merging BAM files for condition {wildcards.condition}",
    shell: """
    samtools merge {output} {input} -@ {threads}
    """


use rule index_bam as index_bam_merged with:
    input:
        bam=OUTPUT_DIR + "/mcc/{condition}/{condition}.bam",
    output:
        bai=OUTPUT_DIR + "/mcc/{condition}/{condition}.bam.bai",
    wildcard_constraints:
        condition="|".join(SAMPLE_GROUPINGS.get_grouping('condition').group_names),
    log: OUTPUT_DIR + "/logs/index_bam_merged/{condition}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/index_bam_merged/{condition}.tsv",
    message: "Indexing merged BAM for condition {wildcards.condition}"


use rule make_bigwigs_mcc_replicates as make_bigwigs_mcc_grouped_condition_norm with:
    input:
        bam=OUTPUT_DIR + "/mcc/{condition}/{condition}.bam",
        bai=OUTPUT_DIR + "/mcc/{condition}/{condition}.bam.bai",
        excluded_regions=OUTPUT_DIR + "/resources/exclusion_regions.bed",
        cis_or_trans_stats=OUTPUT_DIR + "/resources/{condition}_ligation_stats.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/n_cis/{condition}_{viewpoint_group}.bigWig",
    params:
        scale_factor=lambda wc: get_n_cis_scaling_factor(wc),
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    log: OUTPUT_DIR + "/logs/bigwig/{condition}_{viewpoint_group}_n_cis.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwig/{condition}_{viewpoint_group}_n_cis.tsv",
    message: "Generating n_cis normalized bigWig for MCC group {wildcards.condition} and viewpoint group {wildcards.viewpoint_group}",
    container: 'oras://ghcr.io/alsmith151/seqnado_pipeline:latest'
