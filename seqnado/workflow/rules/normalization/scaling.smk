# Endogenous normalization i.e. using total read counts
use rule feature_counts as feature_counts_genome with:
    input:
        bam=expand(OUTPUT_DIR + "/aligned/{sample}.bam", sample=SAMPLE_NAMES),
        bai=expand(OUTPUT_DIR + "/aligned/{sample}.bam.bai", sample=SAMPLE_NAMES),
        annotation=OUTPUT_DIR + "/resources/genomic_bins.saf",
    output:
        counts=OUTPUT_DIR + "/readcounts/feature_counts/read_counts.tsv",
    params:
        options=str(CONFIG.third_party_tools.subread.feature_counts.command_line_arguments.add_include('-F SAF')),
    threads: 
        CONFIG.third_party_tools.subread.feature_counts.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/readcounts/featurecounts/featurecounts.log",
    benchmark: OUTPUT_DIR + "/.benchmark/readcounts/featurecounts/featurecounts.tsv",
    message: "Calculating feature counts for all samples"


# rule setup_for_scaling_factors:
#     input:
#         counts=rules.feature_counts_genome.output.counts,
#     output:
#         formatted_counts=OUTPUT_DIR + "/counts/{group}_formatted_counts.tsv",
#         metadata=OUTPUT_DIR + "/counts/{group}_metadata.tsv",
#     run:
#         counts = format_feature_counts(input[0])
#         counts.to_csv(output[0], sep="\t", index=False)

#         metadata = create_metadata(counts)
#         metadata.to_csv(output[1], sep="\t", index=False, header=False)


# rule calculate_scaling_factors:
#     input:
#         formatted_counts=OUTPUT_DIR + "/counts/{group}_formatted_counts.tsv",
#         metadata=OUTPUT_DIR + "/counts/{group}_metadata.tsv",
#     output:
#         scaling_factors=OUTPUT_DIR + "/resources/{group}_scaling_factors.tsv",
#     container:
#         "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
#     script:
#         "../../scripts/calculate_scaling_factors.R"