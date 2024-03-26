from seqnado.helpers import check_options



rule feature_counts:
    input:
        bam=expand("seqnado_output/aligned/{sample}.bam", sample=SAMPLE_NAMES),
        bai=expand("seqnado_output/aligned/{sample}.bam.bai", sample=SAMPLE_NAMES),
        annotation=config["genome"]["gtf"],
    output:
        counts="seqnado_output/feature_counts/read_counts.tsv",
    params:
        options=check_options(config["featurecounts"]["options"]),
    threads: config["featurecounts"]["threads"]
    resources:
        mem=lambda wildcards, attempt: f"{3 * 2 ** (attempt)}GB",
        runtime="2h",
    log:
        "seqnado_output/logs/readcounts/featurecounts/featurecounts.log",
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
