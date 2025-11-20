from seqnado.helpers import define_time_requested, define_memory_requested

rule feature_counts:
    input:
        bam=expand(OUTPUT_DIR + "/aligned/{sample}.bam", sample=SAMPLE_NAMES),
        bai=expand(OUTPUT_DIR + "/aligned/{sample}.bam.bai", sample=SAMPLE_NAMES),
        annotation=CONFIG.genome.gtf,
    output:
        counts=OUTPUT_DIR + "/readcounts/feature_counts/read_counts.tsv",
    params:
        options=str(CONFIG.third_party_tools.subread.feature_counts.command_line_arguments),
    threads: 
        CONFIG.third_party_tools.subread.feature_counts.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/readcounts/featurecounts/featurecounts.log",
    benchmark: OUTPUT_DIR + "/.benchmarks/readcounts/featurecounts/featurecounts.tsv",
    message: "Running featureCounts to quantify reads for all samples"
    shell: """
    featureCounts \
    -a \
    {input.annotation} \
    -T \
    {threads} \
    --donotsort \
    {params.options} \
    -o \
    {output.counts} \
    {input.bam} \
    > {log} 2>&1
    """
