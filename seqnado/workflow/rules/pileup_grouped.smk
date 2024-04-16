def format_deeptools_options_grouped(wildcards, options):
    is_paired = True

    if not is_paired:
        options = re.sub(r"--extendReads", "", options)
        options = re.sub(r"-e", "", options)
    if not options:
        return ""
    else:
        return options



rule deeptools_make_bigwigs_consensus:
    input:
        bam="seqnado_output/aligned/merged/{sample}.bam",
        bai="seqnado_output/aligned/merged/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/merged/{sample}.bigWig",
    params:
        options=lambda wildcards: format_deeptools_options_grouped(wildcards, config["deeptools"]["bamcoverage"]),
    resources:
        mem="4GB",
        runtime="4h",
    threads:
        config["deeptools"]["threads"]
    log:
        "seqnado_output/logs/bigwigs/{sample}.log",
    shell:
        """
        bamCoverage {params.options} -p {threads} -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """
