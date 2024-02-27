
def get_bam_files_for_merge(wildcards, groups=None):
    return DESIGN.groupby("group").get_group(wildcards.group).fn.tolist()


rule merge_bams:
    input:
        bams=get_bam_files_for_merge,
    output:
        "seqnado_output/consensus_peaks/bam/{group}.bam",
    threads: 8
    log:
        "seqnado_output/consensus_peaks/{group}.log",
    shell:
        """
        samtools merge {output} {input} -@ {threads}
        """


use rule index_bam as index_consensus_bam with:
    input:
        bam="seqnado_output/consensus_peaks/bam/{group}.bam",
    output:
        bai="seqnado_output/consensus_peaks/bam/{group}.bam.bai",
    threads: 8


use rule deeptools_make_bigwigs as deeptools_make_bigwigs_consensus with:
    input:
        bam="seqnado_output/consensus_peaks/bam/{group}.bam",
        bai="seqnado_output/consensus_peaks/bam/{group}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/consensus/{group}.bigWig",
    threads: 8
    log:
        "seqnado_output/logs/consensus_peaks/bigwigs/{group}.log",


rule lanceotron_no_input_consensus:
    input:
        group="seqnado_output/bigwigs/consensus/{group}.bigWig",
    output:
        peaks="seqnado_output/consensus_peaks/{treatment}.bed",

