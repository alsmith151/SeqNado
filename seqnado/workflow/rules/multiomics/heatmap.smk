from seqnado.helpers import define_memory_requested, define_time_requested

from seqnado.outputs.multiomics import get_assay_bigwigs


rule multiomics_heatmap_matrix:
    input:
        rules.gather_bigwigs.output.bw_dir,
    output:
        matrix=OUTPUT_DIR + "multiomics/heatmap/heatmap_matrix.mat.gz",
    params:
        # Use GTF from the first assay config (assuming shared genome)
        gtf=lambda wildcards: LOADED_CONFIGS[ASSAYS[0]]["genome"]["gtf"],
        options=lambda wildcards: str(LOADED_CONFIGS[ASSAYS[0]].get("third_party_tools", {}).get("deeptools", {}).get("compute_matrix", {}).get("command_line_arguments", "")),
        # Collect bigWig files at runtime after assay rules complete
        bigwigs=lambda wildcards: get_assay_bigwigs(wildcards, ASSAYS=ASSAYS, rules=rules),
    threads: lambda wildcards: LOADED_CONFIGS[ASSAYS[0]].get("third_party_tools", {}).get("deeptools", {}).get("compute_matrix", {}).get("threads", 8)
    resources:
        runtime=lambda wildcards, attempt: f"{1 * 2**attempt}h",
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=1),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "multiomics/logs/heatmap/matrix.log",
    benchmark: OUTPUT_DIR + "multiomics/.benchmark/heatmap/matrix.tsv",
    message: "Computing Multiomic heatmap matrix from bigWig files across all assays"
    shell: """
    computeMatrix scale-regions \
    -p {threads} {params.options} \
    --smartLabels \
    --missingDataAsZero \
    -S {params.bigwigs} \
    -R {params.gtf} \
    -o {output.matrix} >> {log} 2>&1
    """


rule multiomics_heatmap_plot:
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
    message: "Generating Multiomic heatmap from matrix"
    shell: """
    plotHeatmap -m {input.matrix} -out {output.heatmap} {params.options}
    """


rule multiomics_heatmap_metaplot:
    input:
        matrix=OUTPUT_DIR + "multiomics/heatmap/heatmap_matrix.mat.gz",
    output:
        metaplot=OUTPUT_DIR + "multiomics/heatmap/metaplot.pdf",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=1),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "multiomics/logs/heatmap/metaplot.log",
    benchmark: OUTPUT_DIR + "multiomics/.benchmark/heatmap/metaplot.tsv",
    message: "Generating Multiomic metaplot from heatmap matrix"
    shell: """
    plotProfile -m {input.matrix} -out {output.metaplot} --perGroup
    """
