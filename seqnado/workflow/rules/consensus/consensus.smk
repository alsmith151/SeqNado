from seqnado.helpers import define_time_requested, define_memory_requested

# Combine bam files for merge

def get_bam_files_for_consensus(wildcards):
    """Get BAM files for merging based on sample names."""
    groups = SAMPLE_GROUPINGS.groupings.get(wildcards.group)
    sample_names = groups.get_samples()
    bam_files = [
        f"seqnado_output/aligned/{sample}.bam" for sample in sample_names
    ]
    return bam_files


rule merge_bams:
    input:
        bams=get_bam_files_for_consensus,
    output:
        temp("seqnado_output/aligned/merged/{group}.bam"),
    threads: CONFIG.third_party_tools.samtools.merge.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: "seqnado_output/logs/merge_bam/{group}.log",
    shell:"""
    samtools merge {output} {input} -@ {threads}
    """


use rule index_bam as index_consensus_bam with:
    input:
        bam="seqnado_output/aligned/merged/{group}.bam",
    output:
        bai=temp("seqnado_output/aligned/merged/{group}.bam.bai"),
    threads: 8



# Pileup for grouped sample

rule deeptools_make_bigwigs_consensus:
    input:
        bam="seqnado_output/aligned/merged/{sample}.bam",
        bai="seqnado_output/aligned/merged/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/merged/{sample}.bigWig",
    params:
        options=lambda wildcards: format_deeptools_options_grouped(wildcards, config["deeptools"]["bamcoverage"]),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads:
        CONFIG.third_party_tools.deeptools.bam_coverage.threads,
    log:
        "seqnado_output/logs/bigwigs/{sample}.log",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        """
        bamCoverage {params.options} -p {threads} -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """

# Peak call for grouped samples

rule lanceotron_no_input_consensus:
    input:
        bigwig="seqnado_output/bigwigs/deeptools/merged/{group}.bigWig",
    output:
        peaks="seqnado_output/peaks/merged/lanceotron/{group}.bed",
        ltron_peaks=temp("seqnado_output/peaks/merged/lanceotron/{group}_L-tron.bed"),
    threads:
        CONFIG.third_party_tools.lanceotron.callpeak.threads
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=10, attempts=attempt, scale=SCALE_RESOURCES),
    params:
        outdir="seqnado_output/peaks/merged/lanceotron",
        options=str(CONFIG.third_party_tools.lanceotron.callpeak.command_line_arguments)
    container:
        "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest"
    log:
        "seqnado_output/logs/lanceotron/{group}.log",
    shell:"""
    lanceotron callPeaks {input.bigwig} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
    cat {output.ltron_peaks} | cut -f 1-3 > {output.peaks}
    """


# Counting for merged samples

rule merged_saf:
    input:
        peaks="seqnado_output/peaks/merged/lanceotron/{group}.bed",
    output:
        saf=temp("seqnado_output/readcounts/featurecounts/{group}.saf"),
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:"""
    awk 'BEGIN{{OFS="\\t"}}{{print $1":"$2"-"$3,$1,$2,$3,"\\*"}}' {input.peaks} > {output.saf}
    """

rule merged_counts:
    input:
        bam=get_bam_files_for_consensus,
        bai=lambda wildcards: [Path(b).with_suffix(".bai") for b in get_bam_files_for_consensus(wildcards)],
        saf=rules.merged_saf.output.saf,
    output:
        counts="seqnado_output/readcounts/featurecounts/{group}_counts.tsv",
    params:
        options=str(CONFIG.third_party_tools.subread.feature_counts.command_line_arguments),
    threads: CONFIG.third_party_tools.subread.feature_counts.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:"seqnado_output/logs/readcounts/featurecounts/{group}_counts.log",
    shell:"""
    featureCounts -a {input.saf} -F SAF -T {threads} --donotsort {params.options} -o {output.counts} {input.bam} > {log} 2>&1
    """

localrules:
    merged_saf