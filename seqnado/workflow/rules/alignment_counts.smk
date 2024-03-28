from seqnado.helpers import check_options

rule feature_counts:
    input:
        bam=expand("seqnado_output/aligned/{sample}.bam", sample=SAMPLE_NAMES),
        bai=expand("seqnado_output/aligned/{sample}.bam.bai", sample=SAMPLE_NAMES),
        annotation=config["genome"]["gtf"],
    output:
        counts="seqnado_output/quantification/feature_counts/read_counts.tsv",
    params:
        options=check_options(config["featurecounts"]["options"]),
    threads: config["featurecounts"]["threads"]
    resources:
        mem=lambda wildcards, attempt: f"{3 * 2 ** (attempt)}GB",
        runtime="2h",
    log:
        "seqnado_output/logs/quantification/featurecounts/featurecounts.log",
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

rule salmon_counts_paired:
    input:
        fq1="seqnado_output/fastq/{sample}_1.fastq.gz", sample=SAMPLE_NAMES,
        fq2="seqnado_output/fastq/{sample}_2.fastq.gz", sample=SAMPLE_NAMES,
    output:
        counts="seqnado_output/quantification/salmon/quant.sf",
    params:
        index=config["salmon"]["index"],
        options=check_options(config["salmon"]["options"]),
    threads: config["salmon"]["threads"]
    resources:
        mem=lambda wildcards, attempt: f"{3 * 2 ** (attempt)}GB",
        runtime="2h",
    log:
        "seqnado_output/logs/readcounts/salmon/salmon.log",
    shell:
        """
        salmon quant -p -t {params.index} {params.options} -1 {input.fq1} -2 {input.fq2} -p {threads} -o seqnado_output/quantification/salmon
        """

rule salmon_counts_single:
    input:
        fq="seqnado_output/fastq/{sample}.fastq.gz", sample=SAMPLE_NAMES,
    output:
        counts="seqnado_output/quantification/salmon/quant.sf",
    params:
        index=config["salmon"]["index"],
        options=check_options(config["salmon"]["options"]),
    threads: config["salmon"]["threads"]
    resources:
        mem=lambda wildcards, attempt: f"{3 * 2 ** (attempt)}GB",
        runtime="2h",
    log:
        "seqnado_output/logs/readcounts/salmon/salmon.log",
    shell:
        """
        salmon quant -t {params.index} {params.options} -r {input.fq} -p {threads} -o seqnado_output/quantification/salmon
        """

ruleorder: feature_counts > salmon_counts_paired > salmon_counts_single
