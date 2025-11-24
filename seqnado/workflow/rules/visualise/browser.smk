
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
        assay=str(CONFIG.assay) if hasattr(CONFIG, 'assay') else None,
        genes=str(CONFIG.assay_config.plotting.genes) if CONFIG.assay_config.plotting and hasattr(CONFIG.assay_config.plotting, 'genes') else None,
        outdir=OUTPUT_DIR + "/genome_browser_plots/",
        regions=str(CONFIG.assay_config.plotting.coordinates) if CONFIG.assay_config.plotting and hasattr(CONFIG.assay_config.plotting, 'coordinates') else None,
        plotting_format=str(CONFIG.assay_config.plotting.file_format) if CONFIG.assay_config.plotting and hasattr(CONFIG.assay_config.plotting, 'file_format') else None,
    resources:
        mem="1.5GB",
         runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "library://asmith151/plotnado/plotnado:latest"
    log: OUTPUT_DIR + "/logs/visualise/plotnado.log",
    benchmark: OUTPUT_DIR + "/.benchmarks/visualise/plotnado.tsv",
    message: "Generating genome browser visualisations with Plotnado"
    script:
        "../../scripts/run_plotnado.py"


        