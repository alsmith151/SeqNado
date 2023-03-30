
rule heatmap_matrix:
    input:
       bigwig=expand(
            "seqnado_output/bigwigs/{method}/{ip}.bigWig",
        ip=[
        f"{row.sample}_{row.antibody}"
                for row in FASTQ_SAMPLES.design.itertuples()
            ],
            method=PILEUP_METHODS,
        ), 
        gtf = config["genome"]["gtf"],
    output:
        matrix = "seqnado_output/heatmap/matrix/matrix.mat.gz",
    params: 
        options = config["heatmap"]["options"],
    threads: 8
    shell: "computeMatrix reference-point --referencePoint TSS -a 3000 -b 3000 -p {threads} --smartLabels --missingDataAsZero {params.options} -S {input.bigwig} -R {input.gtf} -o {output.matrix}"

rule heatmap_plot:
    input:
        matrix = rules.heatmap_matrix.output.matrix,
    output:
        heatmap = "seqnado_output/heatmap/heatmap.pdf",
    params:
        colorMap = config["heatmap"]["colorMap"],
    shell: "plotHeatmap --colorMap {params.colorMap} --boxAroundHeatmaps no -m {input.matrix} -out {output.heatmap}"
