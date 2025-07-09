
rule get_fasta:
    input:
        peaks="seqnado_output/peaks/{method}/{sample}.bed",
    output:
        fasta="seqnado_output/motifs/fasta/{sample}.fa",
        bed=temp("motifs/fasta/{sample}_clean.bed"),
    params:
        genome=config["genome"]["fasta"],
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log: "seqnado_output/logs/motifs/fasta/{sample}.log",
    benchmark: "seqnado_output/benchmarks/motifs/fasta/{sample}.benchmark",
    shell:
        """
    cat {input.peaks} | cut -f 1-3 > {output.bed} &&
    bedtools getfasta -fullHeader -fi {params.genome} -bed {output.bed} -fo {output.fasta}
    """


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
    log: "seqnado_output/logs/motifs/meme/{sample}.log",
    benchmark: "seqnado_output/benchmarks/motifs/meme/{sample}.benchmark",
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
        homer_dir="seqnado_output/motifs/homer/{sample}/",
        homer_params=config["homer"]["homer_params"],
        homer_bg=config["homer"]["homer_bg"],
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log: "seqnado_output/logs/motifs/homer/{sample}.log",
    benchmark: "seqnado_output/benchmarks/motifs/homer/{sample}.benchmark",
    shell:
        """
        findMotifsGenome.pl {input.fasta} {params.homer_bg} {params.homer_dir} {params.homer_params}
        """


ruleorder: motif_meme_chip > motif_homer
