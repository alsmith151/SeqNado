
rule generate_plotnado_visualisation:
    input:
        data=[
            OUTPUT.bigwigs,
            OUTPUT.peaks,
        ],
    output:
        plots=OUTPUT.genome_browser_plots,
        template=OUTPUT_DIR + "/genome_browser_plots/template.toml",
    params:
        assay=ASSAY,
        genes=CONFIG.assay_config.plotting.genes,
        outdir=OUTPUT_DIR + "/genome_browser_plots/",
        regions=CONFIG.assay_config.plotting.regions,
        plotting_format=CONFIG.assay_config.plotting.file_format,
    resources:
        mem="1.5GB",
         runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container:
        "library://asmith151/plotnado/plotnado:latest"
    script:
        "../scripts/run_plotnado.py"


        