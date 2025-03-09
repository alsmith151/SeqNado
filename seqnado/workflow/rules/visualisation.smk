

rule generate_plotnado_visualisation:
    input:
        data=[
            OUTPUT.bigwigs,
            OUTPUT.peaks,
        ],
    output:
        plots=OUTPUT.plots,
        template="seqnado_output/genome_browser_plots/template.toml",
    params:
        assay=ASSAY,
        genes=config['plotting_genes'],
        outdir="seqnado_output/genome_browser_plots/",
        regions=config['plotting_coordinates'],
        plotting_format=config['plotting_format'],
    resources:
        mem="1.5GB",
         runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container:
        "library://asmith151/plotnado/plotnado:latest"
    script:
        "../scripts/run_plotnado.py"


        