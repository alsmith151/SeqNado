rule heatmap_matrix:
    input:
       bigwig="seqnado_output/bigwigs/{method}/{ip}.bigWig",
       gtf = config["genome"]["gtf"],
    output:
        matrix=temp("seqnado_output/heatmap/{method}/matrix/{ip}.mat.gz"),
    params: 
        options = utils.check_options(config["heatmap"]["options"]),
    threads: config["deeptools"]["threads"],
    log: "seqnado_output/logs/heatmap/{method}/matrix/{ip}.log",
    shell: """computeMatrix reference-point --referencePoint TSS -a 3000 -b 3000 \
    -p {threads} --smartLabels --missingDataAsZero {params.options} \
    -S {input.bigwig} \
    -R {input.gtf} -o {output.matrix} >> {log} 2>&1"""

rule heatmap_plot:
    input:
        matrix = rules.heatmap_matrix.output.matrix,
    output:
        heatmap="seqnado_output/heatmap/{method}/heatmap/{ip}.pdf",
    params:
        colormap = config["heatmap"]["colormap"],
    log: "seqnado_output/logs/heatmap/{method}/plot/{ip}.log",
    shell: "plotHeatmap --colorMap {params.colormap} --boxAroundHeatmaps no -m {input.matrix} -out {output.heatmap}"


