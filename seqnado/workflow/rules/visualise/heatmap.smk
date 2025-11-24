from seqnado.helpers import define_memory_requested, define_time_requested


rule heatmap_matrix:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=ScaleMethod.UNSCALED
        ),
    output:
        matrix=temp(OUTPUT_DIR + "/heatmap/heatmap_matrix.mat.gz"),
    params:
        gtf=CONFIG.genome.gtf,
        options=str(CONFIG.third_party_tools.deeptools.compute_matrix.command_line_arguments),
    threads: CONFIG.third_party_tools.deeptools.compute_matrix.threads
    resources:
        runtime=lambda wildcards, attempt: f"{1 * 2**attempt}h",
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/heatmap/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/matrix.tsv",
    message: "Computing heatmap matrix from bigWig files"
    shell: """
    computeMatrix scale-regions -p {threads} {params.options} --smartLabels --missingDataAsZero -S {input.bigwigs} -R {params.gtf} -o {output.matrix} >> {log} 2>&1
    """


rule heatmap_plot:
    input:
        matrix=OUTPUT_DIR + "/heatmap/heatmap_matrix.mat.gz",
    output:
        heatmap=OUTPUT_DIR + "/heatmap/heatmap.pdf",
    params:
        options=str(CONFIG.third_party_tools.deeptools.plot_heatmap.command_line_arguments),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/heatmap/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/heatmap.tsv",
    message: "Generating heatmap from matrix"
    shell: """
    plotHeatmap -m {input.matrix} -out {output.heatmap} {params.options}
    """


rule heatmap_metaplot:
    input:
        matrix=OUTPUT_DIR + "/heatmap/heatmap_matrix.mat.gz",
    output:
        metaplot=OUTPUT_DIR + "/heatmap/metaplot.pdf",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/heatmap/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/metaplot.tsv",
    message: "Generating metaplot from heatmap matrix"
    shell: """
    plotProfile -m {input.matrix} -out {output.metaplot} --perGroup
    """
