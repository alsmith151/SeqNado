from seqnado.helpers import check_options

rule lanceotron_no_input_consensus:
    input:
        bigwig="seqnado_output/bigwigs/deeptools/merged/{group}.bigWig",
    output:
        peaks="seqnado_output/peaks/merged/lanceotron/{group}.bed",
    threads: 8
    resources:
         runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=10, attempts=attempt, scale=SCALE_RESOURCES),,
    params:
        outdir="seqnado_output/peaks/merged/lanceotron",
        options=check_options(config["lanceotron"]["callpeak"])
    container:
        "library://asmith151/seqnado/seqnado_extra:latest"
    log:
        "seqnado_output/logs/lanceotron/{group}.log",
    shell:
        """
        lanceotron callPeaks {input.bigwig} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
        cat {params.outdir}/{wildcards.group}_L-tron.bed | cut -f 1-3 > {output.peaks}
        """