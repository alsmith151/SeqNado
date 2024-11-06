

rule generate_plotnado_visualisation:
    input:
        data=[
            OUTPUT.bigwigs,
            OUTPUT.peaks,
        ],
        regions=config["plotnado"]["regions"],
    output:
        plots=OUTPUT.plots,
        template="seqnado_output/genome_browser_plots/template.toml",
    params:
        assay=ASSAY,
        genome=config["genome"]["name"],
    container:
        "library://asmith151/plotnado/plotnado:latest"
    script:
        "../scripts/run_plotnado.py"


        