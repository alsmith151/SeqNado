import seqnado.utils as utils


def get_bigwigs_for_heatmap():
    bigwigs = []
    if ASSAY == "ChIP":
        for row in FASTQ_SAMPLES.design.itertuples():
            bigwigs.append(
                f"seqnado_output/bigwigs/deeptools/{row.sample}_{row.antibody}.bigWig"
            ),
            sample_heatmap.append(
                f"{row.sample}_{row.antibody}"
            )
    elif ASSAY == "ATAC":
        for row in FASTQ_SAMPLES.design.itertuples():
            bigwigs.append(
                f"seqnado_output/bigwigs/deeptools/{row.sample}.bigWig"
            ),
            sample_heatmap.append(
                f"{row.sample}"
            )
    
    return bigwigs, sample_heatmap



rule heatmap_matrix:
    input:
       bigwig=get_bigwigs_for_heatmap(),
       gtf = config["genome"]["gtf"],
    output:
        matrix=temp("seqnado_output/heatmap/{method}/matrix/{sample_heatmap}.mat.gz"),
    params: 
        options = utils.check_options(config["heatmap"]["options"]),
    threads: config["deeptools"]["threads"],
    log: "seqnado_output/logs/heatmap/{method}/matrix/{sample_heatmap}.log",
    shell: """computeMatrix reference-point --referencePoint TSS -a 3000 -b 3000 \
    -p {threads} --smartLabels --missingDataAsZero {params.options} \
    -S {input.bigwig} \
    -R {input.gtf} -o {output.matrix} >> {log} 2>&1"""

rule heatmap_plot:
    input:
        matrix = rules.heatmap_matrix.output.matrix,
    output:
        heatmap="seqnado_output/heatmap/{method}/heatmap/{sample_heatmap}.pdf",
    params:
        colormap = config["heatmap"]["colormap"],
    log: "seqnado_output/logs/heatmap/{method}/plot/{sample_heatmap}.log",
    shell: "plotHeatmap --colorMap {params.colormap} --boxAroundHeatmaps no -m {input.matrix} -out {output.heatmap}"


