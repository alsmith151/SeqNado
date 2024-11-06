

rule generate_plotnado_visualisation:
    input:
        data=[
            OUTPUT.bigwigs,
            OUTPUT.peaks,
        ],
        regions=config['plotting_coordinates'],
    output:
        plots=OUTPUT.plots,
        template="seqnado_output/genome_browser_plots/template.toml",
    params:
        assay=ASSAY,
        genome=config["genome"]["name"].split('_')[0] if "_" in config["genome"]["name"] else config["genome"]["name"],
    container:
        "library://asmith151/plotnado/plotnado:latest"
    script:
        "../scripts/run_plotnado.py"


        