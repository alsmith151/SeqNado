
rule generate_plotnado_visualisation:
    input:
        data=[
            OUTPUT.bigwig_files,
            OUTPUT.peak_files,
        ],
    output:
        plots=OUTPUT.genome_browser_plots,
        template=OUTPUT_DIR + "/genome_browser_plots/template.toml",
    params:
        assay=CONFIG.assay,
        genes=CONFIG.assay_config.plotting.genes,
        outdir=OUTPUT_DIR + "/genome_browser_plots/",
        regions=CONFIG.assay_config.plotting.coordinates,
        plotting_format=CONFIG.assay_config.plotting.file_format,
    resources:
        mem="1.5GB",
         runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    # container: "library://asmith151/plotnado/plotnado:latest"
    log: OUTPUT_DIR + "/logs/visualise/plotnado.log",
    benchmark: OUTPUT_DIR + "/.benchmarks/visualise/plotnado.tsv",
    message: "Generating genome browser visualisations with Plotnado"
    shell:
        """
        echo "Running plotnado for genome browser visualisations" >> {log} 2>&1
        echo "params: {params}" >> {log} 2>&1
        touch {output.template}
        touch {output.plots}
        echo "Completed plotnado visualisations" >> {log} 2>&1
        """
    # script:
    #     "../../scripts/run_plotnado.py"


        