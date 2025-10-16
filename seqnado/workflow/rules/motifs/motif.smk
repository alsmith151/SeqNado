rule motif_meme_chip:
    input:
        fasta="seqnado_output/motifs/fasta/{sample}.fa",
    output:
        meme="seqnado_output/motifs/meme/{sample}/meme-chip.html",
    params:
        meme_dir="seqnado_output/motifs/meme/{sample}/",
        meme_chip_params=config["meme"]["meme_chip_params"],
        meme_chip_db=config["meme"]["meme_chip_db"],
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/motifs/meme/{sample}.log",
    shell:
        """
        meme-chip -oc {params.meme_dir} {params.meme_chip_db} {params.meme_chip_params} {input.fasta}
        """


rule motif_homer:
    input:
        fasta="seqnado_output/motifs/fasta/{sample}.fa",
    output:
        homer="seqnado_output/motifs/homer/{sample}/homerResults.html",
    params:
        genome=CONFIG.genome.fasta,
        homer_dir="seqnado_output/motifs/homer/{sample}/",
        homer_params=CONFIG.third_party_tools.homer.find_motifs_genome.command_line_arguments,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/motifs/homer/{sample}.log",
    shell:
        """
        findMotifsGenome.pl {input.fasta} {params.genome} {params.homer_dir} {params.homer_params}
        """


ruleorder: motif_meme_chip > motif_homer
