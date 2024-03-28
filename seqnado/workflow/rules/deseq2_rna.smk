rule deseq2_report_rnaseq:
    input:
        counts="seqnado_output/readcounts/feature_counts/read_counts.tsv",
        qmd=f"deseq2_{PROJECT_NAME}.qmd".replace(" ", ""),
    output:
        deseq2=f"deseq2_{PROJECT_NAME}.html".replace(" ", ""),
        size_factors="seqnado_output/resources/all_normalisation_factors.json"
    log:
        "seqnado_output/logs/deseq2/deseq2.log",
    container:
        "library://asmith151/seqnado/seqnado_report:latest"
    shell:
        """
        input_file=$(realpath "{input.qmd}")
        base_dir=$(dirname $input_file)
        cd "$base_dir"
        quarto render {input.qmd} --no-cache --output {output.deseq2} --log {log}
        """


localrules:
    deseq2_report_rnaseq,
