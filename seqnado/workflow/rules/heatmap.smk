import seqnado.utils

if ASSAY == "ChIP":
    prefix = SAMPLE_NAMES_IP 
elif ASSAY == "RNA":
    prefix = [x + y for x in SAMPLE_NAMES for y in ["_plus", "_minus"]]
else:
    prefix = SAMPLE_NAMES


rule heatmap_matrix:
    input:
        bigwigs=expand("seqnado_output/bigwigs/{method}/{sample}.bigWig", sample=prefix, method=config["pileup_method"]),
    output:
        matrix=temp("seqnado_output/heatmap/heatmap_matrix.mat.gz"),
    params: 
        gtf = config["genome"]["gtf"],
        options = utils.check_options(config["heatmap"]["options"]),
    threads: 
        config["deeptools"]["threads"],
    resources:
        time="02:00:00",
        mem_mb=lambda wildcards, attempt: 16000 * 2**attempt,
    log: 
        "seqnado_output/logs/heatmap/matrix.log",
    shell: """computeMatrix scale-regions --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --regionBodyLength 5000 --smartLabels --missingDataAsZero -p {threads} {params.options} -S {input.bigwigs} -R {params.gtf} -o {output.matrix} >> {log} 2>&1"""

rule heatmap_plot:
    input:
        matrix="seqnado_output/heatmap/heatmap_matrix.mat.gz",
    output:
        heatmap="seqnado_output/heatmap/heatmap.pdf",
    params:
        colormap = utils.check_options(config["heatmap"]["colormap"]),
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * 2**attempt,
    log: 
        "seqnado_output/logs/heatmap/heatmap.log",
    shell: """plotHeatmap -m {input.matrix} -out {output.heatmap} --colorMap {params.colormap} --boxAroundHeatmaps no"""


rule heatmap_metaplot:
    input:
        matrix="seqnado_output/heatmap/heatmap_matrix.mat.gz",
    output:
        metaplot="seqnado_output/heatmap/metaplot.pdf",
    params:
        colormap = utils.check_options(config["heatmap"]["colormap"]),
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * 2**attempt,
    log: 
        "seqnado_output/logs/heatmap/metaplot.log",
    shell: """plotProfile -m {input.matrix} -out {output.metaplot} --perGroup"""
