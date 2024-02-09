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
        bigwig="seqnado_output/consensus_peaks/{sample}.bigWig",
    threads: 8


rule lanceotron_no_input_consensus:
    input:
        group="seqnado_output/consensus_peaks/{group}.bigWig",
    output:
        peaks="seqnado_output/consensus_peaks/{group}.bed",
    log:
        "seqnado_output/consensus_peaks/{group}.log",
    params:
        options=seqnado.utils.check_options(config["lanceotron"]["callpeak"]),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
    threads: 1
    container:
        "library://asmith151/seqnado/seqnado_extra:latest"
    resources:
        mem_mb=10_1000,
        time="0-06:00:00",
    shell:
        """
        lanceotron callPeaks {input.group} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
        cat {params.outdir}/{wildcards.group}_L-tron.bed | cut -f 1-3 > {output.peaks}
        """


use rule bed_to_bigbed as bed_to_bigbed_consensus with:
    input:
        bed="seqnado_output/consensus_peaks/{sample}.bed",
    output:
        bigbed="seqnado_output/consensus_peaks/{sample}.bigBed",
    threads: 1
