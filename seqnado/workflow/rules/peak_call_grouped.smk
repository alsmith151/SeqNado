from seqnado.helpers import check_options

rule lanceotron_no_input_consensus:
    input:
        bigwig="seqnado_output/bigwigs/deeptools/merged/{group}.bigWig",
    output:
        peaks="seqnado_output/peaks/merged/lanceotron/{group}.bed",
    threads: 8
    resources:
        runtime=lambda wildcards, attempt: f"{6 * 2 ** (attempt - 1)}h",
        mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt - 1)}GB",
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