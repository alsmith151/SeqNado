rule extract_ligation_stats:
    input:
        bam="seqnado_output/mcc/replicates/{sample}/{sample}.bam",
    output:
        stats="seqnado_output/resources/{sample}_ligation_stats.json"
    container: 'oras://ghcr.io/alsmith151/seqnado_pipeline:latest'
    shell:
        """
        mccnado extract-ligation-stats {input.bam} {output.stats} 
        """




use rule extract_ligation_stats as extract_ligation_stats_merged with:
    input:
        bam="seqnado_output/mcc/{group}/{group}.bam",
    output:
        stats="seqnado_output/resources/{group}_ligation_stats.json",
    log:
        "seqnado_output/logs/extract_ligation_stats_merged/{group}.log",