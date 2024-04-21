import pathlib


rule deseq2_report_rnaseq:
    input:
        counts="seqnado_output/readcounts/feature_counts/read_counts.tsv",
        qmd=f"deseq2_{PROJECT_NAME}.qmd".replace(" ", ""),
        yml="seqnado_output/resources/deseq2_params.yml"
    output:
        deseq2=f"deseq2_{PROJECT_NAME}.html".replace(" ", ""),
    log:
        "seqnado_output/logs/deseq2/deseq2.log",
    container:
        "library://asmith151/seqnado/seqnado_report:latest"
    shell:
        """
        input_file=$(realpath "{input.qmd}")
        base_dir=$(dirname $input_file)
        cd "$base_dir"
        quarto render {input.qmd} --no-cache --output {output.deseq2} --log {log} --execute-params {input.yml}
        """



rule deseq2_params:
    output:
        yml="seqnado_output/resources/deseq2_params.yml"
    params:
        spikein_genes=["AmpR_seq", "Cas9_5p_seq", "Cas9_3p_seq"],
        size_factors_out="seqnado_output/resources/all_normalisation_factors.json",
        de_dir=str(pathlib.Path(rules.deseq2_report_rnaseq.output.deseq2).parent),
        counts=rules.deseq2_report_rnaseq.input.counts,
    container: None
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
