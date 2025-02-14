from seqnado.helpers import check_options, define_time_requested, define_memory_requested

rule create_consensus_saf:
    input:
        peaks="seqnado_output/peaks/merged/lanceotron/{group}.bed",
    output:
        consensus_saf="seqnado_output/readcounts/{group}.saf",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/readcounts/{group}_saf.log",
    shell:"""
        awk 'BEGIN{OFS="\t"}{print $1":"$2"-"$3,$1,$2,$3,"+"}' {input.peaks} > {output.consensus_saf} > {log} 2>&1
        """


rule feature_counts_consensus:
    input:
        peaks="seqnado_output/peaks/merged/lanceotron/{group}.bed",
        bam=expand("seqnado_output/aligned/{sample}.bam", sample=SAMPLE_NAMES),
        bai=expand("seqnado_output/aligned/{sample}.bam.bai", sample=SAMPLE_NAMES),
    output:
        saf=temp("seqnado_output/readcounts/{group}.saf"),
        counts="seqnado_output/readcounts/feature_counts/{group}_counts.tsv",
    params:
        options="--primary --ignoreDup --maxMOp 20"
    threads: config["featurecounts"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/readcounts/featurecounts/featurecounts.log",
    shell:
        """
        awk 'BEGIN{{OFS="\t"}}{{print $1":"$2"-"$3,$1,$2,$3,"+"}}' {input.peaks} > {output.saf} &&
        featureCounts -a {output.saf} -T {threads} {params.options} -o {output.counts} {input.bam} > {log} 2>&1
        """
