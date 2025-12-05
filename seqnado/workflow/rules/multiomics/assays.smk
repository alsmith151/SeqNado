
rule run_assay:
    """
    Execute a single assay pipeline with outputs in seqnado_output/{assay}/
    
    This runs the main Snakefile with modified output paths.
    Runs as a regular rule (not checkpoint) to enable parallel execution of multiple assays.
    """
    output:
        complete = OUTPUT_DIR + "{assay}/logs/.complete",
        report = OUTPUT_DIR + "{assay}/seqnado_report.html"
    params:
        assay = lambda wildcards: wildcards.assay,
        config_file = lambda wildcards: str(Path(ASSAY_CONFIGS[wildcards.assay]['path']).resolve()),
        metadata = lambda wildcards: str(Path(f"metadata_{wildcards.assay}.csv").resolve()),
        main_snakefile = str(Path(workflow.basedir) / "Snakefile"),
        output_dir = OUTPUT_DIR + "{assay}",
        project_root = str(Path.cwd()),
        cores_for_assay = CORES_PER_ASSAY
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 259200  # 72 hours
    log:
        OUTPUT_DIR + "{assay}/logs/seqnado.log"
    message:
        "Running {params.assay} assay pipeline with config {params.config_file}"
    shell:
        """
        # Create output directory
        mkdir -p {params.output_dir}/logs
        
        echo "================================================================" | tee {log}
        echo "Starting {params.assay} pipeline at $(date)" | tee -a {log}
        echo "Config: {params.config_file}" | tee -a {log}
        echo "Metadata: {params.metadata}" | tee -a {log}
        echo "Output: {params.output_dir}" | tee -a {log}
        echo "================================================================" | tee -a {log}
        
        # Run the main SeqNado Snakefile for this assay
        # Create assay-specific .snakemake directory to prevent locking conflicts
        mkdir -p .snakemake/assay_{params.assay}

        # Use --directory to isolate .snakemake metadata in assay-specific subdirectory
        # The .snakemake dir will be created at .snakemake/assay_{{params.assay}}/.snakemake
        cd .snakemake/assay_{params.assay} && \
        snakemake \
            --snakefile {params.main_snakefile} \
            --configfile {params.config_file} \
            --config output_dir=../../{params.output_dir} metadata={params.metadata} \
            --cores {params.cores_for_assay} \
            {WORKFLOW_ARGS} \
            2>&1 | tee -a ../../{log}

        # Mark completion (from .snakemake/assay_{{params.assay}}/ so need ../../ prefix)
        touch ../../{output.complete}
        
        echo "================================================================" | tee -a ../../{log}
        echo "Completed {params.assay} pipeline at $(date)" | tee -a ../../{log}
        echo "Report: ../../{output.report}" | tee -a ../../{log}
        echo "================================================================" | tee -a ../../{log}
        """