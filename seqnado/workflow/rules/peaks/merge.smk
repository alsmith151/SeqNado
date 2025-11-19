
rule lanceotron_no_input_consensus:
    input:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/{group}.bigWig",
    output:
        peaks=OUTPUT_DIR + "/peaks/merged/lanceotron/{group}.bed",
        ltron_peaks=temp(OUTPUT_DIR + "/peaks/merged/lanceotron/{group}_L-tron.bed"),
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
    log:
        OUTPUT_DIR + "/logs/lanceotron/{group}.log",
    shell:"""
    lanceotron callPeaks {input.bigwig} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
    cat {output.ltron_peaks} | cut -f 1-3 > {output.peaks}
    """