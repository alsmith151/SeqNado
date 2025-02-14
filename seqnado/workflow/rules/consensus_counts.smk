from seqnado.helpers import check_options, define_time_requested, define_memory_requested


def get_bam_files_for_counts(wildcards):
    from seqnado.design import NormGroups
    norm_groups = NormGroups.from_design(DESIGN, subset_column="merge")
    sample_names = norm_groups.get_grouped_samples(wildcards.group)
    print(f"Sample names: {sample_names}")
    bam_files = [
        f"seqnado_output/aligned/{sample}.bam" for sample in sample_names
    ]
    return bam_files


def get_bai_files_for_counts(wildcards):
    from seqnado.design import NormGroups
    norm_groups = NormGroups.from_design(DESIGN, subset_column="merge")

    sample_names = norm_groups.get_grouped_samples(wildcards.group)
    bai_files = [
        f"seqnado_output/aligned/{sample}.bam.bai" for sample in sample_names
    ]
    return bai_files


rule feature_counts_consensus:
    input:
        peaks="seqnado_output/peaks/merged/lanceotron/{group}.bed",
        bam=get_bam_files_for_counts,
        bai=get_bai_files_for_counts,
    output:
        saf=temp("seqnado_output/readcounts/{group}.saf"),
        counts="seqnado_output/readcounts/feature_counts_{group}/read_counts.tsv",
    params:
        options=config["featurecounts"]["options"],
    threads: config["featurecounts"]["threads"],
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log:"seqnado_output/logs/readcounts/featurecounts/featurecounts_{group}.log",
    shell:
        """
        awk 'BEGIN{{OFS="\t"}}{{print $1":"$2"-"$3,$1,$2,$3,"+"}}' {input.peaks} > {output.saf} &&
        featureCounts -a {output.saf} -T {threads} {params.options} -o {output.counts} {input.bam} > {log} 2>&1
        """
