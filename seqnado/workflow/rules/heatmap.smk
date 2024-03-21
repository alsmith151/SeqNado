from seqnado.helpers import check_options, get_scale_method

if ASSAY == "ChIP":
    prefix = SAMPLE_NAMES_IP
elif ASSAY == "RNA":
    prefix = [x + y for x in SAMPLE_NAMES for y in ["_plus", "_minus"]]
else:
    prefix = SAMPLE_NAMES


rule heatmap_matrix:
    input:
        bigwigs=expand(
            "seqnado_output/bigwigs/deeptools/{method}/{sample}.bigWig",
            method=get_scale_method(config),
            sample=prefix,
        ),
    output:
        matrix=temp("seqnado_output/heatmap/heatmap_matrix.mat.gz"),
    params:
        gtf=config["genome"]["gtf"],
        options=check_options(config["heatmap"]["options"]),
    threads: config["deeptools"]["threads"]
    resources:
        runtime=lambda wildcards, attempt: f"{1 * 2**attempt}h",
        mem=lambda wildcards, attempt: f"{4 * 2**attempt}GB",
    log:
        "seqnado_output/logs/heatmap/matrix.log",
    shell:
        """computeMatrix scale-regions -p {threads} {params.options} --smartLabels --missingDataAsZero -S {input.bigwigs} -R {params.gtf} -o {output.matrix} >> {log} 2>&1"""


rule heatmap_plot:
    input:
        matrix="seqnado_output/heatmap/heatmap_matrix.mat.gz",
    output:
        heatmap="seqnado_output/heatmap/heatmap.pdf",
    params:
        colormap=check_options(config["heatmap"]["colormap"]),
    resources:
        mem=lambda wildcards, attempt: f"{2 * 2**attempt}GB",
    log:
        "seqnado_output/logs/heatmap/heatmap.log",
    shell:
        """plotHeatmap -m {input.matrix} -out {output.heatmap} --colorMap {params.colormap} --boxAroundHeatmaps no"""


rule heatmap_metaplot:
    input:
        matrix="seqnado_output/heatmap/heatmap_matrix.mat.gz",
    output:
        metaplot="seqnado_output/heatmap/metaplot.pdf",
    resources:
        mem=lambda wildcards, attempt: f"{2 * 2**attempt}GB"
    log:
        "seqnado_output/logs/heatmap/metaplot.log",
    shell:
        """plotProfile -m {input.matrix} -out {output.metaplot} --perGroup"""
