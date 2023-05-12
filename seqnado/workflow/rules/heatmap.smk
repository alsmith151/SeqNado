from seqnado import utils

rule heatmap_matrix:
    input:
        bigwig="seqnado_output/bigwigs/{wildcards.method}/{wildcards.sample}.bigWig",
        gtf = config["genome"]["gtf"],
    output:
        matrix=temp("seqnado_output/heatmap/{wildcards.method}/matrix/{wildcards.sample}.mat.gz"),
    params: 
        options = utils.check_options(config["heatmap"]["options"]),
    threads: 
        config["deeptools"]["threads"],
    resources:
        time="02:00:00",
    log: 
        "seqnado_output/logs/heatmap/{method}/{sample}.log",
    shell: 
        """computeMatrix \
           reference-point \
           --referencePoint TSS \
           -a 3000 -b 3000 \
           -p {threads} \
           --smartLabels \
           --missingDataAsZero \
           {params.options} \
           -S {input.bigwig} \
           -R {input.gtf} \
           -o {output.matrix} >> {log} 2>&1"""

rule heatmap_plot:
    input:
        matrix = "seqnado_output/heatmap/{wildcards.method}/matrix/{wildcards.sample}.mat.gz",
    output:
        heatmap="seqnado_output/heatmap/{method}/{sample}.png",
    params:
        colormap = utils.check_options(config["heatmap"]["colormap"]),
    resources:
        mem_mb=1024 * 2,
    log: 
        "seqnado_output/logs/heatmap/{method}/plot/{sample}.log",
    shell: 
        """plotHeatmap -m {input.matrix} \
        -out {output.heatmap} \
        --colorMap {params.colormap} \
        --boxAroundHeatmaps no"""
