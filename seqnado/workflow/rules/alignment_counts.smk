from seqnado.helpers import define_time_requested, define_memory_requested

rule feature_counts:
    input:
        bam=expand("seqnado_output/aligned/{sample}.bam", sample=SAMPLE_NAMES),
        bai=expand("seqnado_output/aligned/{sample}.bam.bai", sample=SAMPLE_NAMES),
        annotation=CONFIG.genome.gtf,
    output:
        counts="seqnado_output/readcounts/feature_counts/read_counts.tsv",
    params:
        options=str(CONFIG.third_party_tools.subread.feature_counts.command_line_arguments),
    threads: 
        CONFIG.third_party_tools.subread.feature_counts.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/readcounts/featurecounts/featurecounts.log",
    shell:
        """
        featureCounts \
        -a \
        {input.annotation} \
        -T \
        {threads} \
        --donotsort \
        {params.options} \
        -o \
        {output.counts} \
        {input.bam} \
        > {log} 2>&1
        """

rule salmon_counts_paired:
    input:
        fq1="seqnado_output/fastqs/{sample}_1.fastq.gz",
        fq2="seqnado_output/fastqs/{sample}_2.fastq.gz",
    output:
        counts="seqnado_output/readcounts/salmon/salmon_{sample}/quant.sf",
        out_dir=temp(directory("seqnado_output/readcounts/salmon/salmon_{sample}"))
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.salmon.quant.command_line_arguments),
    threads:
        CONFIG.third_party_tools.salmon.quant.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/readcounts/salmon/salmon_{sample}.log",
    shell:
        """
        salmon quant -i {params.index} {params.options} -1 {input.fq1} -2 {input.fq2} -p {threads} -o {output.out_dir}
        """

rule salmon_counts_single:
    input:
        fq="seqnado_output/fastqs/{sample}.fastq.gz"
    output:
        counts="seqnado_output/readcounts/salmon/salmon_{sample}/quant.sf",
        out_dir=temp(directory("seqnado_output/readcounts/salmon/salmon_{sample}"))
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.salmon.quant.command_line_arguments),
    threads: CONFIG.third_party_tools.salmon.quant.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/readcounts/salmon/salmon_{sample}.log",
    shell:
        """
        salmon quant -i {params.index} {params.options} -r {input.fq} -p {threads} -o {output.out_dir}
        """

rule get_salmon_counts:
    input:
        count_dirs=expand("seqnado_output/readcounts/salmon/salmon_{sample}", sample=SAMPLE_NAMES),
        counts=expand("seqnado_output/readcounts/salmon/salmon_{sample}/quant.sf", sample=SAMPLE_NAMES)
    output:
        count_table="seqnado_output/readcounts/salmon/salmon_counts.csv"
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/readcounts/salmon/salmon_counts.log"
    script:
        "../scripts/get_salmon_counts.py"

ruleorder: feature_counts > salmon_counts_paired > salmon_counts_single
