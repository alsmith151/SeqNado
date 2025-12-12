from seqnado.helpers import define_memory_requested, define_time_requested


def get_all_assay_bigwigs(wildcards):
    """Collect bigwig files from all individual assays."""
    import glob
    bigwigs = []
    for assay in ASSAYS:
        # Construct the path to bigwigs from each assay's output directory
        assay_output_dir = f"{OUTPUT_DIR}{assay}"
        # Pattern matches deeptools unscaled bigwigs (same as single-assay heatmap)
        pattern = f"{assay_output_dir}/bigwigs/deeptools/unscaled/*.bigWig"
        assay_bws = sorted(glob.glob(pattern))
        if assay_bws:
            print(f"[DEBUG] Found {len(assay_bws)} bigwigs for {assay}: {assay_bws}")
        else:
            print(f"[DEBUG] No bigwigs found for {assay} at pattern: {pattern}")
        bigwigs.extend(assay_bws)

    if not bigwigs:
        raise ValueError(
            f"No bigwig files found for any assay. Checked patterns:\n" +
            "\n".join([f"  - {OUTPUT_DIR}{assay}/bigwigs/deeptools/unscaled/*.bigWig" for assay in ASSAYS])
        )

    return bigwigs


rule multiassay_heatmap_matrix:
    input:
        bigwigs=get_all_assay_bigwigs,
        summary=OUTPUT_DIR + "multi_assay_summary.txt"
    output:
        matrix=OUTPUT_DIR + "multiomics/heatmap/heatmap_matrix.mat.gz",
    params:
        # Use GTF from the first assay config (assuming shared genome)
        gtf=lambda wildcards: LOADED_CONFIGS[ASSAYS[0]]["genome"]["gtf"],
        options=lambda wildcards: str(LOADED_CONFIGS[ASSAYS[0]].get("third_party_tools", {}).get("deeptools", {}).get("compute_matrix", {}).get("command_line_arguments", "")),
    threads: lambda wildcards: LOADED_CONFIGS[ASSAYS[0]].get("third_party_tools", {}).get("deeptools", {}).get("compute_matrix", {}).get("threads", 8)
    resources:
        runtime=lambda wildcards, attempt: f"{1 * 2**attempt}h",
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=1),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "multiomics/logs/heatmap/matrix.log",
    benchmark: OUTPUT_DIR + "multiomics/.benchmark/heatmap/matrix.tsv",
    message: "Computing multi-assay heatmap matrix from bigWig files across all assays"
    shell: """
    computeMatrix scale-regions -p {threads} {params.options} --smartLabels --missingDataAsZero -S {input.bigwigs} -R {params.gtf} -o {output.matrix} >> {log} 2>&1
    """


rule multiassay_heatmap_plot:
    input:
        matrix=OUTPUT_DIR + "multiomics/heatmap/heatmap_matrix.mat.gz",
    output:
        heatmap=OUTPUT_DIR + "multiomics/heatmap/heatmap.pdf",
    params:
        options=lambda wildcards: str(LOADED_CONFIGS[ASSAYS[0]].get("third_party_tools", {}).get("deeptools", {}).get("plot_heatmap", {}).get("command_line_arguments", "")),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=1),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "multiomics/logs/heatmap/heatmap.log",
    benchmark: OUTPUT_DIR + "multiomics/.benchmark/heatmap/heatmap.tsv",
    message: "Generating multi-assay heatmap from matrix"
    shell: """
    plotHeatmap -m {input.matrix} -out {output.heatmap} {params.options}
    """


rule multiassay_heatmap_metaplot:
    input:
        matrix=OUTPUT_DIR + "multiomics/heatmap/heatmap_matrix.mat.gz",
    output:
        metaplot=OUTPUT_DIR + "multiomics/heatmap/metaplot.pdf",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=1),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "multiomics/logs/heatmap/metaplot.log",
    benchmark: OUTPUT_DIR + "multiomics/.benchmark/heatmap/metaplot.tsv",
    message: "Generating multi-assay metaplot from heatmap matrix"
    shell: """
    plotProfile -m {input.matrix} -out {output.metaplot} --perGroup
    """