def get_bam_files_for_merge(wildcards):
    # Find all sample names that belong to the group
    gb = DESIGN.to_dataframe().groupby("merge")
    samples = gb.groups[wildcards.group]

    # Get the file names for the samples
    return expand("seqnado_output/aligned/{sample}.bam", sample=samples)


rule merge_bams:
    input:
        bams=get_bam_files_for_merge,
    output:
        "seqnado_output/consensus_peaks/{group}.bam",
    threads: 8
    log:
        "seqnado_output/consensus_peaks/{group}.log",
    shell:
        """
        samtools merge {output} {input} -@ {threads}
        """


use rule index_bam as index_consensus_bam with:
    input:
        bam="seqnado_output/consensus_peaks/{sample}.bam",
    output:
        bai="seqnado_output/consensus_peaks/{sample}.bam.bai",
    threads: 8


use rule deeptools_make_bigwigs as deeptools_make_bigwigs_consensus with:
    input:
        bam="seqnado_output/consensus_peaks/{sample}.bam",
        bai="seqnado_output/consensus_peaks/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/consensus_peaks/{sample}.bw",
    threads: 8


use rule lanceotron_no_input as lanceotron_no_input_consensus with:
    input:
        treatment="seqnado_output/consensus_peaks/{treatment}.bw",
    output:
        peaks="seqnado_output/consensus_peaks/{treatment}.bed",
    log:
        "seqnado_output/consensus_peaks/{treatment}.log",
