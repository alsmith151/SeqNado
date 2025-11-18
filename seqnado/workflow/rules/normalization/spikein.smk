
if CONFIG.assay_config.spikein.method == SpikeInMethod.ORLANDO:

    rule calculate_normalisation_factors:
        input:
            lambda wc: expand(
                rules.split_bam.output.stats,
                sample=SAMPLE_NAMES_IP + SAMPLE_NAMES_CONTROL,
            ),
        output:
            normalisation_table="seqnado_output/resources/{group}_normalisation_factors.tsv",
            normalisation_factors="seqnado_output/resources/{group}_normalisation_factors.json",
        log:
            "seqnado_output/logs/normalisation_factors_{group}.log",
        script:
            "../scripts/calculate_spikein_norm_orlando.py"

elif CONFIG.assay_config.spikein.method == SpikeInMethod.WITH_INPUT:

    rule calculate_normalisation_factors:
        input:
            lambda wc: expand(
                rules.split_bam.output.stats,
                sample=SAMPLE_NAMES_IP + SAMPLE_NAMES_CONTROL,
            ),
            design="seqnado_output/metadata.csv",
        output:
            normalisation_table="seqnado_output/resources/{group}_normalisation_factors.tsv",
            normalisation_factors="seqnado_output/resources/{group}_normalisation_factors.json",
        log:
            "seqnado_output/logs/normalisation_factors_{group}.log",
        script:
            "../scripts/calculate_spikein_norm_factors.py"