from seqnado.helpers import define_time_requested, define_memory_requested


# Pileup for grouped sample

rule deeptools_make_bigwigs_consensus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/{group}.bigWig",
    params:
       options=lambda wildcards: format_deeptools_options(wildcards, str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments)),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads:
        CONFIG.third_party_tools.deeptools.bam_coverage.threads,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwigs/deeptools/merged/{group}.tsv",
    message: "Making bigWig with deeptools for merged sample {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        bamCoverage {params.options} -p {threads} -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """