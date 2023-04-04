import seqnado.utils as utils


def get_bigwigs_for_heatmap():
    bigwigs = []
    if ASSAY == "ChIP":
        for row in FASTQ_SAMPLES.design.itertuples():
            bigwigs.append(
                f"seqnado_output/bigwigs/deeptools/{row.sample}_{row.antibody}.bigWig"
            )
    elif ASSAY == "ATAC":
        for row in FASTQ_SAMPLES.design.itertuples():
            bigwigs.append(f"seqnado_output/bigwigs/deeptools/{row.sample}.bigWig")
    
    return bigwigs



rule heatmap_matrix:
    input:
       bigwig=get_bigwigs_for_heatmap(),
       gtf = config["genome"]["gtf"],
    output:
        matrix = "seqnado_output/heatmap/matrix/matrix.mat.gz",
    params: 
        options = utils.check_options(config["heatmap"]["options"]),
    threads: 8
    shell: "computeMatrix reference-point --referencePoint TSS -a 3000 -b 3000 -p {threads} --smartLabels --missingDataAsZero {params.options} -S {input.bigwig} -R {input.gtf} -o {output.matrix}"

rule heatmap_plot:
    input:
        matrix = rules.heatmap_matrix.output.matrix,
    output:
        heatmap = "seqnado_output/heatmap/heatmap.pdf",
    params:
        colormap = config["heatmap"]["colormap"],
    shell: "plotHeatmap --colorMap {params.colormap} --boxAroundHeatmaps no -m {input.matrix} -out {output.heatmap}"
