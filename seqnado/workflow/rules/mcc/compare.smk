rule make_aggregate_bigwigs:
    input:
        bigwigs=expand(OUTPUT_DIR + "/bigwigs/mcc/n_cis/{group}_{viewpoint_group}.bigWig",
               group=SAMPLE_GROUPINGS.get_grouping('condition').group_names,
               viewpoint_group=VIEWPOINT_TO_GROUPED_VIEWPOINT.values()),
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/aggregated-using-mean/{group}_{viewpoint_group}.bigWig",
    params:
        options=str(CONFIG.third_party_tools.bamnado.bigwig_aggregate.command_line_arguments),
    log: OUTPUT_DIR + "/logs/bigwig/{group}_{viewpoint_group}_aggregated-using-mean.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwig/{group}_{viewpoint_group}_aggregated-using-mean.tsv",
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    message: "Generating consensus bigWig for MCC group {wildcards.group} and viewpoint group {wildcards.viewpoint_group}",
    shell: """
    bamnado \
    bigwig-aggregate \
    --bigwigs {input.bigwigs} \
    -o {output.bigwig} \
    -m mean \
    {params.options} \
    > {log} 2>&1
    """


rule make_comparison_bigwigs:
    input:
        bw1=OUTPUT_DIR + "/bigwigs/mcc/aggregated-using-mean/{group1}_{viewpoint_group}.bigWig",
        bw2=OUTPUT_DIR + "/bigwigs/mcc/aggregated-using-mean/{group2}_{viewpoint_group}.bigWig",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/subtractions/{group1}_vs_{group2}_{viewpoint_group}.bigWig",
    params:
        options=str(CONFIG.third_party_tools.bamnado.bigwig_compare.command_line_arguments),
    log: OUTPUT_DIR + "/logs/bigwig/{group1}_vs_{group2}_{viewpoint_group}_subtraction.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwig/{group1}_vs_{group2}_{viewpoint_group}_subtraction.tsv",
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    message: "Generating comparison bigWig for MCC groups {wildcards.group1} vs {wildcards.group2} and viewpoint group {wildcards.viewpoint_group}",
    shell: """
    bamnado \
    bigwig-compare \
    --bw1 {input.bw1} \
    --bw2 {input.bw2} \
    -o {output.bigwig} \
    -c subtraction \
    {params.options} \
    > {log} 2>&1
    """


