
# Define input and output files so that seqnado isnt required in the rule

def get_assay_bigwigs(wildcards):
    """Get all bigwigs from assay-specific 'all' rules."""
    bigwigs = []
    peak_files = []
    for assay in ASSAYS:
        rule_name = f"{assay}_all"
        inputs = getattr(rules, rule_name).input

        # Handle both single files and collections of files
        if isinstance(inputs, str):
            if inputs.endswith(".bigWig"):
                bigwigs.append(inputs)
        else:
            # InputFiles is iterable
            for file in inputs:
                if isinstance(file, str) and file.endswith(".bigWig"):
                    bigwigs.append(file)

        if isinstance(inputs, str):
            if inputs.endswith(".bed") 
                peak_files.append(inputs)
        else:
            # InputFiles is iterable
            for file in inputs:
                if isinstance(file, str) and file.endswith(".bed"):
                    peak_files.append(file)

        input_files= bigwigs + peak_files
    return input_files

plot_files=OUTPUT.genome_browser_plots

rule generate_plotnado_visualisation:
    input:
        data=get_assay_bigwigs,
    output:
        plots=plot_files,
        template=OUTPUT_DIR + "multiomics/genome_browser_plots/template.toml",
    params:
        assay=str(CONFIG.assay) if hasattr(CONFIG, 'assay') else None,
        genes=str(CONFIG.genome.genes) if CONFIG.assay_config.plot_with_plotnado and hasattr(CONFIG.genome, 'genes') else None,
        outdir=OUTPUT_DIR + "multiomics/genome_browser_plots/",
        regions=str(CONFIG.assay_config.plotting.coordinates) if CONFIG.assay_config.plotting and hasattr(CONFIG.assay_config.plotting, 'coordinates') else None,
        plotting_format=str(CONFIG.assay_config.plotting.file_format) if CONFIG.assay_config.plotting and hasattr(CONFIG.assay_config.plotting, 'file_format') else None,
    resources:
        mem="1.5GB",
         runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "library://asmith151/plotnado/plotnado:latest"
    log: OUTPUT_DIR + "multiomics/logs/visualise/plotnado.log",
    benchmark: OUTPUT_DIR + "multiomics/.benchmark/visualise/plotnado.tsv",
    message: "Generating genome browser visualisations with Plotnado"
    script:
        "../../scripts/run_plotnado.py"
