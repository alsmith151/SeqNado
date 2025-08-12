from pathlib import Path
from seqnado.helpers import define_time_requested, define_memory_requested

rule merged_saf:
    input:
        peaks="seqnado_output/peaks/merged/lanceotron/{group}.bed",
    output:
        saf=temp("seqnado_output/readcounts/featurecounts/{group}.saf"),
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
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
    log:"seqnado_output/logs/readcounts/featurecounts/{group}_counts.log",
    shell:"""
    featureCounts -a {input.saf} -F SAF -T {threads} --donotsort {params.options} -o {output.counts} {input.bam} > {log} 2>&1
    """

localrules:
    merged_saf