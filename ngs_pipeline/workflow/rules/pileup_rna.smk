import re


def split_options_and_scale_factor():

    options = [op.strip() for op in config["deeptools"]["bamcoverage"].split()]

    if "--scaleFactor" in options:
        option_index = options.index("--scaleFactor")
        scale_factor = options.pop(option_index + 1)
        options.pop(option_index)

    else:
        scale_factor = 1

    return {"options": " ".join(options), "scale_factor": scale_factor}


rule deeptools_make_bigwigs_rna_plus:
    input:
        bam="aligned_and_filtered/{sample}.bam",
        bai="aligned_and_filtered/{sample}.bam.bai",
        filtering="flags/{sample}.filtering.complete.sentinel",
    output:
        bigwig="bigwigs/deeptools/{sample}_plus.bigWig",
    params:
        **split_options_and_scale_factor(),
    threads: config["deeptools"]["threads"]
    log:
        "logs/pileups/deeptools/{sample}.log",
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --filterRNAstrand forward -p {threads} --scaleFactor {params.scale_factor} > {log} 2>&1
        """


rule deeptools_make_bigwigs_rna_minus:
    input:
        bam="aligned_and_filtered/{sample}.bam",
        bai="aligned_and_filtered/{sample}.bam.bai",
        filtering="flags/{sample}.filtering.complete.sentinel",
    output:
        bigwig="bigwigs/deeptools/{sample}_minus.bigWig",
    params:
        **split_options_and_scale_factor(),
    threads: config["deeptools"]["threads"]
    log:
        "logs/pileups/deeptools/{sample}.log",
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --filterRNAstrand reverse -p {threads} --scaleFactor -{params.scale_factor} > {log} 2>&1
        """
