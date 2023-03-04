import seqnado.utils as utils

rule feature_counts:
    input:
        bam=expand("aligned/{sample}.bam", sample=SAMPLE_NAMES),
        bai=expand("aligned/{sample}.bam.bai", sample=SAMPLE_NAMES),
        annotation=config["genome"]["annotation"],
    output:
        counts=f"feature_counts/read_counts.tsv",
    params:
        options=utils.check_options(config["featurecounts"]["options"]),
    threads: config["featurecounts"]["threads"]
    log:
        "logs/readcounts/featurecounts/featurecounts.log",
    shell:
        """
        featureCounts \
        -a \
        {input.annotation} \
        -T \
        {threads} \
        {params.options} \
        -o \
        {output.counts} \
        {input.bam} \
        > {log} 2>&1
        """
