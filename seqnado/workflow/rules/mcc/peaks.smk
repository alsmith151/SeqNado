

rule call_mcc_peaks: # TODO: ensure that we're using the GPU queue
    input:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/unscaled/{group}_{viewpoint_group}.bigWig",
    output:
        peaks=OUTPUT_DIR + "/peaks/lanceotron-mcc/{group}_{viewpoint_group}.bed",
    log:
        OUTPUT_DIR + "/logs/call_mcc_peaks/{group}_{viewpoint_group}.log",
    params:
        options=str(CONFIG.third_party_tools.lanceotron_mcc.call_peaks.command_line_arguments),
    container: None
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        gpu=1,
    shell:
        """
        apptainer exec \
        --nv \
        library://asmith151/lanceotron/lanceotron-mcc:latest \
        lanceotron-mcc \
        call-mcc-peaks \
        --bigwig {input.bigwig} \
        --outfile {output.peaks} \
        {params.options} > {log} 2>&1
        """

