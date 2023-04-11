
localrules: deseq2_report_rnaseq

rule deseq2_report_rnaseq:
    input:
        counts="seqnado_output/feature_counts/read_counts.tsv",
        qmd=f"DESeq2_{config['deseq2']['project_id']}.qmd",
    output:
        deseq2=f"DESeq2_{config['deseq2']['project_id']}.html",
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