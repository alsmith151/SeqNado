rule extract_ligation_stats:
    input:
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam",
    output:
        stats=OUTPUT_DIR + "/resources/{sample}_ligation_stats.json"
    container: 'oras://ghcr.io/alsmith151/seqnado_pipeline:latest'
    shell:
        """
        mccnado extract-ligation-stats {input.bam} {output.stats} 
        """




use rule extract_ligation_stats as extract_ligation_stats_merged with:
    input:
        bam=OUTPUT_DIR + "/mcc/{group}/{group}.bam",
    output:
        stats=OUTPUT_DIR + "/resources/{group}_ligation_stats.json",
    log:
        OUTPUT_DIR + "/logs/extract_ligation_stats_merged/{group}.log",