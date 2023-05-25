from seqnado import utils

def get_norm_factor(wildcards):
    with open("seqnado_output/normalisation_factors.json") as f:
        norm_factors = json.load(f)
    return norm_factors[wildcards.sample]

def format_deeptools_bamcoverage_options(wildcards):
    import re

    norm = get_norm_factor(wildcards)

    options = utils.check_options(config["deeptools"]["bamcoverage"])


    if "--normalizeUsing" in options:
        options = re.sub("--normalizeUsing [a-zA-Z]+", "", options)
    
    if "--scaleFactor" in options:
        options = re.sub("--scaleFactor [0-9\.]+", "", options)

    options += f" --scaleFactor {norm}"

    return options

def format_homer_make_bigwigs_options(wildcards):
    import re

    norm = int(get_norm_factor(wildcards) * 1e7)

    options = utils.check_options(config["homer"]["makebigwig"])

    if "-norm" in options:
        options = re.sub("-scale [0-9.]+", "", options)

    options += f" -norm {norm}"

    return options


rule calculate_normalisation_factors:
    input:
        bam_ref = expand("seqnado_output/aligned/{sample}.bam", sample=SAMPLE_NAMES),
        bam_spikein = expand("seqnado_output/aligned/spikein/{sample}.bam", sample=SAMPLE_NAMES),
        bam_spikein_index = expand("seqnado_output/aligned/spikein/{sample}.bam.bai", sample=SAMPLE_NAMES),
        design = rules.save_design.output[0],
    output:
        normalisation_table = "seqnado_output/normalisation_factors.tsv",
        normalisation_factors = "seqnado_output/normalisation_factors.json"
    container:
        None
    log:
        "seqnado_output/logs/normalisation_factors.log"
    script:
        "../scripts/calculate_spikein_norm_factors.py"

use rule deeptools_make_bigwigs as deeptools_make_bigwigs_norm with:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
        normalisation_factors="seqnado_output/normalisation_factors.json",
    params:
        options =lambda wildcards: format_deeptools_bamcoverage_options(wildcards)

use rule homer_make_bigwigs as homer_make_bigwigs_norm with:
    input:
        homer_tag_directory="seqnado_output/tag_dirs/{sample}",
        normalisation_factors="seqnado_output/normalisation_factors.json",
    params:
        genome_name=config["genome"]["name"],
        genome_chrom_sizes=config["genome"]["chromosome_sizes"],
        options=lambda wc: format_homer_make_bigwigs_options(wc),
        outdir="seqnado_output/bigwigs/homer/",
        temp_bw=lambda wc, output: output.homer_bigwig.replace(".bigWig", ".ucsc.bigWig")


    