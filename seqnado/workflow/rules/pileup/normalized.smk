def get_norm_factor_spikein(wildcards, negative=False):
    import json

    norm_file = OUTPUT_DIR + "/resources/normalisation_factors.json"
    with open(norm_file, 'r') as f:
        norm_factors = json.load(f)

    scale_factor = float(norm_factors.get(wildcards.sample, 1.0))

    if negative:
        scale_factor = -scale_factor

    return scale_factor

def get_scaling_factor(wildcards, scaling_file):
    import pandas as pd

    df = pd.read_csv(scaling_file, sep="\t")
    scale = df.loc[df['sample'] == wildcards.sample, 'scale_factor'].values[0]

    return float(scale)
    
rule deeptools_make_bigwigs_scale:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + "/resources/{get_group_for_sample(wc , INPUT_FILES)}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/csaw/{sample}.bigWig",
    params:
        scale=lambda wc: get_scaling_factor(
            wc,
            OUTPUT_DIR + "/resources/{get_group_for_sample(wc , INPUT_FILES)}_scaling_factors.tsv",
        ),
        options=lambda wc: format_deeptools_bamcoverage_options(wc)
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),    
    log: OUTPUT_DIR + "/logs/pileups/deeptools/scaled/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/deeptools/scaled/{sample}.tsv",
    message: "Making scaled bigWig with deeptools for sample {wildcards.sample}"
    shell: """
    bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options} > {log} 2>&1
    """


use rule deeptools_make_bigwigs_scale as deeptools_make_bigwigs_spikein with:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + "/resources/{get_group_for_sample(wc , INPUT_FILES)}_normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/spikein/{sample}.bigWig",
    params:
        options=lambda wildcards: format_deeptools_bamcoverage_options(wildcards),
        scale=get_norm_factor_spikein,
    log: OUTPUT_DIR + "/logs/pileups/deeptools/spikein/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/deeptools/spikein/{sample}.tsv",
    message: "Making spike-in normalized bigWig with deeptools for sample {wildcards.sample}"

rule deeptools_make_bigwigs_rna_spikein_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + "/resources/{get_group_for_sample(wc , INPUT_FILES)}_normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/spikein/{sample}_plus.bigWig",
    params:
        options=lambda wildcards: format_deeptools_bamcoverage_options(wildcards),
        scale=get_norm_factor_spikein,
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/pileups/deeptools/spikein/{sample}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/deeptools/spikein/{sample}_plus.tsv",
    message: "Making plus strand spike-in normalized bigWig with deeptools for sample {wildcards.sample}"
    shell: """
    bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} --scaleFactor {params.scale} {params.options} --filterRNAstrand forward > {log} 2>&1
    """


rule deeptools_make_bigwigs_rna_spikein_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + "/resources/{get_group_for_sample(wc , INPUT_FILES)}_normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/spikein/{sample}_minus.bigWig",
    params:
        options=lambda wildcards: format_deeptools_bamcoverage_options(wildcards),
        scale=lambda wc: get_norm_factor_spikein(wc, negative=True),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/pileups/deeptools/spikein/{sample}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/deeptools/spikein/{sample}_minus.tsv",
    message: "Making minus strand spike-in normalized bigWig with deeptools for sample {wildcards.sample}"
    shell: """
    bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} --scaleFactor {params.scale} {params.options} --filterRNAstrand reverse > {log} 2>&1
    """