from seqnado.workflow.helpers.mcc import get_n_cis_scaling_factor, get_mcc_bam_files_for_merge


rule make_bigwigs_mcc_replicates:
    input:
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam",
        bai=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam.bai",
        excluded_regions=OUTPUT_DIR + "/resources/exclusion_regions.bed",
        cis_or_trans_stats=OUTPUT_DIR + "/resources/replicates/{sample}_ligation_stats.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/replicates/{sample}_{viewpoint_group}.bigWig"
    params:
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
        scale_factor=lambda wc: get_n_cis_scaling_factor(wc, OUTPUT_DIR=OUTPUT_DIR),
    log: OUTPUT_DIR + "/logs/bigwig/{sample}_{viewpoint_group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwig/{sample}_{viewpoint_group}.tsv",
    message: "Generating bigWig for MCC sample {wildcards.sample} and viewpoint group {wildcards.viewpoint_group}",
    shell: """
    bamnado \
    bam-coverage \
    -b {input.bam} \
    -o {output.bigwig} \
    --scale-factor {params.scale_factor} \
    --blacklisted-locations {input.excluded_regions} \
    --min-mapq 0 \
    --read-group {wildcards.viewpoint_group} \
    {params.options} > {log} 2>&1
    """
        


rule merge_mcc_bams:
    input:
        bams=lambda wc: get_mcc_bam_files_for_merge(wc, SAMPLE_GROUPINGS=SAMPLE_GROUPINGS, OUTPUT_DIR=OUTPUT_DIR),
    output:
        OUTPUT_DIR + "/mcc/{group}/{group}.bam",
    threads: CONFIG.third_party_tools.samtools.merge.threads
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/merge_bam/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/merge_bam/{group}.tsv",
    message: "Merging BAM files for group {wildcards.group}",
    shell: """
    samtools merge {output} {input} -@ {threads}
    """


use rule index_bam as index_bam_merged with:
    input:
        bam=OUTPUT_DIR + "/mcc/{group}/{group}.bam",
    output:
        bai=OUTPUT_DIR + "/mcc/{group}/{group}.bam.bai",
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    log: OUTPUT_DIR + "/logs/index_bam_merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/index_bam_merged/{group}.tsv",
    message: "Indexing merged BAM for group {wildcards.group}"


use rule make_bigwigs_mcc_replicates as make_bigwigs_mcc_grouped_norm with:
    input:
        bam=OUTPUT_DIR + "/mcc/{group}/{group}.bam",
        bai=OUTPUT_DIR + "/mcc/{group}/{group}.bam.bai",
        excluded_regions=OUTPUT_DIR + "/resources/exclusion_regions.bed",
        cis_or_trans_stats=OUTPUT_DIR + "/resources/{group}_ligation_stats.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/n_cis/{group}_{viewpoint_group}.bigWig",
    params:
        scale_factor=lambda wc: get_n_cis_scaling_factor(wc, OUTPUT_DIR=OUTPUT_DIR),
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    log: OUTPUT_DIR + "/logs/bigwig/{group}_{viewpoint_group}_n_cis.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwig/{group}_{viewpoint_group}_n_cis.tsv",
    message: "Generating n_cis normalized bigWig for MCC group {wildcards.group} and viewpoint group {wildcards.viewpoint_group}",
    container: 'oras://ghcr.io/alsmith151/seqnado_pipeline:latest'


use rule make_bigwigs_mcc_replicates as make_bigwigs_mcc_grouped_raw with:
    input:
        bam=OUTPUT_DIR + "/mcc/{group}/{group}.bam",
        bai=OUTPUT_DIR + "/mcc/{group}/{group}.bam.bai",
        excluded_regions=OUTPUT_DIR + "/resources/exclusion_regions.bed",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/unscaled/{group}_{viewpoint_group}.bigWig"
    params:
        scale_factor=1,
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    log: OUTPUT_DIR + "/logs/bigwig/{group}_{viewpoint_group}_unscaled.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwig/{group}_{viewpoint_group}_unscaled.tsv",
    message: "Generating unscaled bigWig for MCC group {wildcards.group} and viewpoint group {wildcards.viewpoint_group}",


rule confirm_mcc_bigwigs_generated:
    input:
        expand(OUTPUT_DIR + "/bigwigs/mcc/replicates/{sample}_{viewpoint_group}.bigWig", sample=SAMPLE_NAMES, viewpoint_group=VIEWPOINT_TO_GROUPED_VIEWPOINT.values()),
        expand(OUTPUT_DIR + "/bigwigs/mcc/n_cis/{group}_{viewpoint_group}.bigWig", group=SAMPLE_GROUPINGS.get_grouping('consensus').group_names, viewpoint_group=VIEWPOINT_TO_GROUPED_VIEWPOINT.values()),
        expand(OUTPUT_DIR + "/bigwigs/mcc/unscaled/{group}_{viewpoint_group}.bigWig", group=SAMPLE_GROUPINGS.get_grouping('consensus').group_names, viewpoint_group=VIEWPOINT_TO_GROUPED_VIEWPOINT.values()),
    output:
        touch(OUTPUT_DIR + "/bigwigs/mcc/.mcc_bigwigs_generated.txt"),
    message: "Confirming all MCC bigWigs have been generated"
    
