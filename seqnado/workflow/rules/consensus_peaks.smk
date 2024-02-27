def get_bam_files_for_merge(wildcards, groups=None):
    return DESIGN.groupby("group").get_group(wildcards.group).fn.tolist()


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
