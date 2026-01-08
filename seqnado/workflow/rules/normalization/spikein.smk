from seqnado import SpikeInMethod

if CONFIG.assay_config.spikein and CONFIG.assay_config.spikein.method == SpikeInMethod.ORLANDO:

    rule calculate_normalisation_factors:
        input:
            lambda wc: expand(
                rules.split_bam.output.stats,
                sample=SAMPLE_NAMES,
            ),
        output:
            normalisation_table=OUTPUT_DIR + "/resources/orlando/normalisation_factors.tsv",
            normalisation_factors=OUTPUT_DIR + "/resources/orlando/normalisation_factors.json",
        log: OUTPUT_DIR + "/logs/normalisation_factors.log",
        benchmark: OUTPUT_DIR + "/.benchmark/normalisation_factors.tsv",
        message: "Calculating normalisation factors"
        script:
            "../../scripts/calculate_spikein_norm_orlando.py"

elif CONFIG.assay_config.spikein and CONFIG.assay_config.spikein.method == SpikeInMethod.WITH_INPUT:

    rule calculate_normalisation_factors:
        input:
            lambda wc: expand(
                rules.split_bam.output.stats,
                sample=SAMPLE_NAMES,
            ),
            design=OUTPUT_DIR + "/metadata.csv",
        output:
            normalisation_table=OUTPUT_DIR + "/resources/with_input/normalisation_factors.tsv",
            normalisation_factors=OUTPUT_DIR + "/resources/with_input/normalisation_factors.json",
        log: OUTPUT_DIR + "/logs/normalisation_factors.log",
        benchmark: OUTPUT_DIR + "/.benchmark/normalisation_factors.tsv",
        message: "Calculating normalisation factors"
        script:
            "../../scripts/calculate_spikein_norm_factors.py"

elif CONFIG.assay_config.spikein and CONFIG.assay_config.spikein.method == SpikeInMethod.DESEQ2:

    rule calculate_normalisation_factors:
        input:
            counts=OUTPUT_DIR + "/readcounts/feature_counts/read_counts.tsv",
            design=OUTPUT_DIR + "/metadata.csv",
        output:
            normalisation_table=OUTPUT_DIR + "/resources/deseq2/normalisation_factors.tsv",
            normalisation_factors=OUTPUT_DIR + "/resources/deseq2/normalisation_factors.json",
        params:
            spikein_genes=CONFIG.assay_config.spikein.control_genes or [],
            bam_dir=OUTPUT_DIR + "/aligned",
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        log: OUTPUT_DIR + "/logs/normalisation_factors.log",
        benchmark: OUTPUT_DIR + "/.benchmark/normalisation_factors.tsv",
        message: "Calculating DESeq2 spike-in normalisation factors"
        script:
            "../../scripts/calculate_spikein_norm_factors_deseq.R"

elif CONFIG.assay_config.spikein and CONFIG.assay_config.spikein.method == SpikeInMethod.EDGER:

    rule calculate_normalisation_factors:
        input:
            counts=OUTPUT_DIR + "/readcounts/feature_counts/read_counts.tsv",
            design=OUTPUT_DIR + "/metadata.csv",
        output:
            normalisation_table=OUTPUT_DIR + "/resources/edger/normalisation_factors.tsv",
            normalisation_factors=OUTPUT_DIR + "/resources/edger/normalisation_factors.json",
        params:
            spikein_genes=CONFIG.assay_config.spikein.control_genes or [],
            bam_dir=OUTPUT_DIR + "/aligned",
            outdir=OUTPUT_DIR + "/resources/edger",
            outfile="normalisation_factors.json",
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        log: OUTPUT_DIR + "/logs/normalisation_factors.log",
        benchmark: OUTPUT_DIR + "/.benchmark/normalisation_factors.tsv",
        message: "Calculating edgeR spike-in normalisation factors"
        script:
            "../../scripts/calculate_spikein_norm_factors_edger.R"
