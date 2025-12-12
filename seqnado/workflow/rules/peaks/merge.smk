
rule lanceotron_no_input_consensus:
    input:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/{group}.bigWig",
    output:
        peaks=OUTPUT_DIR + "/peaks/merged/lanceotron/merged/{group}.bed",
        ltron_peaks=temp(OUTPUT_DIR + "/peaks/merged/lanceotron/merged/{group}_L-tron.bed"),
    threads:
        CONFIG.third_party_tools.lanceotron.call_peaks.threads
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=10, attempts=attempt, scale=SCALE_RESOURCES),
    params:
        outdir=OUTPUT_DIR + "/peaks/merged/lanceotron",
        options=str(CONFIG.third_party_tools.lanceotron.call_peaks.command_line_arguments)
    container:
        "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest"
    log: OUTPUT_DIR + "/logs/lanceotron/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/lanceotron/{group}.tsv",
    message: "Calling peaks with LanceOtron for merged group {wildcards.group}"
    shell: """
    lanceotron callPeaks {input.bigwig} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
    cat {output.ltron_peaks} | cut -f 1-3 > {output.peaks}
    """

use rule macs2_no_input as macs2_no_input_consensus with:
    input:
        treatment = OUTPUT_DIR + "/aligned/merged/{group}.bam",
    output:
        peaks = OUTPUT_DIR + "/peaks/macs2/merged/{group}.bed",
    log: OUTPUT_DIR + "/logs/macs2/merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/macs2/merged/{group}.tsv",
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    message: "Calling peaks with MACS2 for merged group {wildcards.group}"

