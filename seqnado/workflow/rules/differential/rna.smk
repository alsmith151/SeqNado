from pathlib import Path


rule deseq2_report_rnaseq:
    input:
        counts=OUTPUT_DIR + "/readcounts/feature_counts/read_counts.tsv",
        qmd=f"deseq2_{PROJECT_NAME}.qmd".replace(" ", ""),
        yml=OUTPUT_DIR + "/resources/deseq2_params.yml"
    output:
        deseq2=f"deseq2_{PROJECT_NAME}.html".replace(" ", ""),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/deseq2/deseq2.log",
    benchmark: OUTPUT_DIR + "/.benchmark/deseq2/deseq2.tsv",
    message: "Generating DESeq2 report for RNA-seq data",
    shell: """
    input_file=$(realpath "{input.qmd}")
    base_dir=$(dirname $input_file)
    cd "$base_dir"
    quarto render {input.qmd} --no-cache --output {output.deseq2} --log {log} --execute-params {input.yml}
    """


rule deseq2_params:
    output:
        yml=OUTPUT_DIR + "/resources/deseq2_params.yml"
    params:
        spikein_genes=["AmpR_seq", "Cas9_5p_seq", "Cas9_3p_seq"],
        size_factors_out=OUTPUT_DIR + "/resources/all_normalisation_factors.json",
        de_dir=str(Path(rules.deseq2_report_rnaseq.output.deseq2).parent),
        counts=rules.deseq2_report_rnaseq.input.counts,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/deseq2/deseq2_params.log",
    benchmark: OUTPUT_DIR + "/.benchmark/deseq2/deseq2_params.tsv",
    message: "Preparing DESeq2 parameters for RNA-seq analysis",
    run:
        import yaml

        with open(output.yml, "w") as f:
            yaml.dump(
                {
                    "spikein_genes": params.spikein_genes,
                    "size_factors_out": params.size_factors_out,
                    "de_dir": params.de_dir,
                    "counts": params.counts,
                },
                f,
            )



localrules:
    deseq2_report_rnaseq,
    deseq2_params
