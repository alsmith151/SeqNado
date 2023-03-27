
rule heatmap_matrix:
    input:
        bigwig = expand(
            "seqnado_output/bigwigs/{method}/{sample}.bigWig",
            sample=SAMPLE_NAMES,
            method=PILEUP_METHODS,
        ),
        gtf = index=config["genome"]["gtf"],
    output:
        matrix = "seqnado_output/heatmap/matrix/matrix_{sample}.mat.gz",
    threads: 8
    shell: "computeMatrix reference-point --referencePoint TSS -a 3000 -b 3000 -p {threads} --smartLabels --missingDataAsZero -S {input.bigwig} -R {input.gtf} -o {output.matrix}"

rule heatmap_plot:
    input:
        matrix = rules.heatmap_matrix.output.matrix,
    output:
        heatmap = "seqnado_output/heatmap/heatmap/Heatmap_{sample}.png",
        SortedRegions = "seqnado_output/heatmap/SortedRegions/SortedRegions_{sample}.bed",
    params:
        colorMap = config["heatmap"]["colorMap"],
    shell: "plotHeatmap --colorMap {params.colorMap} --boxAroundHeatmaps no -m {input.matrix} -out {output.heatmap} --outFileSortedRegions {output.SortedRegions}"
