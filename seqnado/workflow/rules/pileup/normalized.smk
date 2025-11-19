
rule deeptools_make_bigwigs_scale:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + "/resources/{get_group_for_sample(wc , DESIGN)}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/csaw/{sample}.bigWig",
    params:
        scale=lambda wc: get_scaling_factor(
            wc,
            OUTPUT_DIR + "/resources/{get_group_for_sample(wc , DESIGN)}_scaling_factors.tsv",
        ),
        options=lambda wc: format_deeptools_bamcoverage_options(wc)
    threads: config["deeptools"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),    
    log:
        OUTPUT_DIR + "/logs/pileups/deeptools/scaled/{sample}.log",
    shell:
        "bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options} > {log} 2>&1"


use rule deeptools_make_bigwigs_scale as deeptools_make_bigwigs_spikein with:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + "/resources/{get_group_for_sample(wc , DESIGN)}_normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/spikein/{sample}.bigWig",
    params:
        options=lambda wildcards: format_deeptools_bamcoverage_options(wildcards),
        scale=get_norm_factor_spikein,


rule deeptools_make_bigwigs_rna_spikein_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + "/resources/{get_group_for_sample(wc , DESIGN)}_normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/spikein/{sample}_plus.bigWig",
    params:
        options=lambda wildcards: format_deeptools_bamcoverage_options(wildcards),
        scale=get_norm_factor_spikein,
    threads: config["deeptools"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        OUTPUT_DIR + "/logs/pileups/deeptools/spikein/{sample}_plus.log",
    shell:
        "bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} --scaleFactor {params.scale} {params.options} --filterRNAstrand forward > {log} 2>&1"


rule deeptools_make_bigwigs_rna_spikein_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + "/resources/{get_group_for_sample(wc , DESIGN)}_normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/spikein/{sample}_minus.bigWig",
    params:
        options=lambda wildcards: format_deeptools_bamcoverage_options(wildcards),
        scale=lambda wc: get_norm_factor_spikein(wc, negative=True),
    threads: config["deeptools"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        OUTPUT_DIR + "/logs/pileups/deeptools/spikein/{sample}_minus.log",
    shell:
        "bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} --scaleFactor {params.scale} {params.options} --filterRNAstrand reverse > {log} 2>&1"