

def get_n_cis_scaling_factor(wc):
    import json 
    from pathlib import Path

    # Can either extract the stats for the sample or the group
    if hasattr(wc, "group"):
        stats_file = OUTPUT_DIR + "/resources/{wc.group}_ligation_stats.json"
    else:
        # If not, then use the sample
        stats_file = OUTPUT_DIR + "/resources/{wc.sample}_ligation_stats.json"
    
    # Create Path object and ensure the file exists
    stats_path = Path(stats_file)
    if not stats_path.exists():
        # raise FileNotFoundError(f"Stats file not found: {stats_file}")
        return 1
        
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
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam",
        bai=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam.bai",
        excluded_regions=OUTPUT_DIR + "/resources/exclusion_regions.bed",
        cis_or_trans_stats=OUTPUT_DIR + "/resources/{sample}_ligation_stats.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/replicates/{sample}_{viewpoint_group}.bigWig"
    params:
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
        scale_factor=lambda wc: get_n_cis_scaling_factor(wc),
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
        
def get_mcc_bam_files_for_merge(wildcards):
    """Get BAM files for merging based on sample names."""
    groups = SAMPLE_GROUPINGS.groupings.get(wildcards.group)
    sample_names = groups.get_samples() if groups else []
    bam_files = [
        OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam" for sample in sample_names
    ]
    return bam_files


rule merge_mcc_bams:
    input:
        bams=get_mcc_bam_files_for_merge,
    output:
        OUTPUT_DIR + "/mcc/{group}/{group}.bam",
    threads: CONFIG.third_party_tools.samtools.merge.threads
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
        scale_factor=lambda wc: get_n_cis_scaling_factor(wc),
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
        bin_size=config['bamnado'].get("bin_size", 10),
        scale_factor=1,
        options=check_options(config["bamnado"]["bamcoverage"]),
    log: OUTPUT_DIR + "/logs/bigwig/{group}_{viewpoint_group}_unscaled.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwig/{group}_{viewpoint_group}_unscaled.tsv",
    message: "Generating unscaled bigWig for MCC group {wildcards.group} and viewpoint group {wildcards.viewpoint_group}",
    


