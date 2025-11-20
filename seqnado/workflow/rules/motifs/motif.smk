rule motif_meme_chip:
    input:
        fasta=OUTPUT_DIR + "/motifs/fasta/{sample}.fa",
    output:
        meme=OUTPUT_DIR + "/motifs/meme/{sample}/meme-chip.html",
    params:
        meme_dir=OUTPUT_DIR + "/motifs/meme/{sample}/",
        meme_chip_params=config["meme"]["meme_chip_params"],
        meme_chip_db=config["meme"]["meme_chip_db"],
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/motifs/meme/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/motifs/meme/{sample}.tsv",
    message: "Running MEME-ChIP motif analysis for sample {wildcards.sample}"
    shell: """
    meme-chip -oc {params.meme_dir} {params.meme_chip_db} {params.meme_chip_params} {input.fasta}
    """


rule motif_homer:
    input:
        fasta=OUTPUT_DIR + "/motifs/fasta/{sample}.fa",
    output:
        homer=OUTPUT_DIR + "/motifs/homer/{sample}/homerResults.html",
    params:
        genome=CONFIG.genome.fasta,
        homer_dir=OUTPUT_DIR + "/motifs/homer/{sample}/",
        homer_params=CONFIG.third_party_tools.homer.find_motifs_genome.command_line_arguments,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/motifs/homer/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/motifs/homer/{sample}.tsv",
    message: "Running HOMER motif analysis for sample {wildcards.sample}"
    shell: """
    findMotifsGenome.pl {input.fasta} {params.genome} {params.homer_dir} {params.homer_params}
    """


ruleorder: motif_meme_chip > motif_homer
