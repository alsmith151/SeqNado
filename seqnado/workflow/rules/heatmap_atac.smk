rule heatmap_matrix:
    input:
       bigwig="seqnado_output/bigwigs/{method}/{sample}.bigWig",
       gtf = config["genome"]["gtf"],
    output:
        matrix=temp("seqnado_output/heatmap/{method}/matrix/{sample}.mat.gz"),
    params: 
        options = utils.check_options(config["heatmap"]["options"]),
    threads: config["deeptools"]["threads"],
    log: "seqnado_output/logs/heatmap/{method}/{sample}.log",
    shell: """computeMatrix reference-point --referencePoint TSS -a 3000 -b 3000 \
    -p {threads} --smartLabels --missingDataAsZero {params.options} \
    -S {input.bigwig} \
    -R {input.gtf} -o {output.matrix} >> {log} 2>&1"""

rule heatmap_plot:
    input:
        matrix = rules.heatmap_matrix.output.matrix,
    output:
        heatmap="seqnado_output/heatmap/{method}/{sample}.png",
    params:
        colormap = config["heatmap"]["colormap"],
    shell: "plotHeatmap --colorMap {params.colormap} --boxAroundHeatmaps no -m {input.matrix} -out {output.heatmap}"
    resources:
        mem_mb=1024 * 10,
    log: "seqnado_output/logs/heatmap/{method}/plot/{sample}.log",
    shell: """plotHeatmap -m {input.matrix} \
    -out {output.heatmap} \
    --colorMap {params.colormap} \
    --boxAroundHeatmaps no"""