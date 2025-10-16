from seqnado.helpers import define_time_requested, define_memory_requested


# Pileup for grouped sample

rule deeptools_make_bigwigs_consensus:
    input:
        bam="seqnado_output/aligned/merged/{sample}.bam",
        bai="seqnado_output/aligned/merged/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/merged/{sample}.bigWig",
    params:
        options=lambda wildcards: format_deeptools_options_grouped(wildcards, config["deeptools"]["bamcoverage"]),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads:
        CONFIG.third_party_tools.deeptools.bam_coverage.threads,
    log:
        "seqnado_output/logs/bigwigs/{sample}.log",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        """
        bamCoverage {params.options} -p {threads} -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """