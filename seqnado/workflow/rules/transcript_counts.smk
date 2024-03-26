from seqnado.helpers import check_options


rule salmon_counts_single:
    input:
        fq="seqnado_output/fastq/{sample}.fastq.gz", sample=SAMPLE_NAMES,
    output:
        counts="seqnado_output/salmon_counts/quant.sf",
    params:
        index=config["salmon"]["index"],
        options=check_options(config["salmon"]["options"]),
    threads: config["salmon"]["threads"]
    resources:
        mem=lambda wildcards, attempt: f"{3 * 2 ** (attempt)}GB",
        runtime="2h",
    log:
        "seqnado_output/logs/readcounts/salmon/salmon.log",
    shell:
        """
        salmon quant -i {params.index} -l A -r {input.fq} -p {threads} -o seqnado_output/salmon_counts
        """
