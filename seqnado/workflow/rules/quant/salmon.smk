from seqnado.helpers import define_time_requested, define_memory_requested

rule salmon_counts_paired:
    input:
        fq1=OUTPUT_DIR + "/fastqs/{sample}_1.fastq.gz",
        fq2=OUTPUT_DIR + "/fastqs/{sample}_2.fastq.gz",
    output:
        counts=OUTPUT_DIR + "/readcounts/salmon/salmon_{sample}/quant.sf",
        out_dir=temp(directory(OUTPUT_DIR + "/readcounts/salmon/salmon_{sample}"))
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
        OUTPUT_DIR + "/logs/readcounts/salmon/salmon_{sample}.log",
    shell:
        """
        salmon quant -i {params.index} {params.options} -1 {input.fq1} -2 {input.fq2} -p {threads} -o {output.out_dir}
        """

rule salmon_counts_single:
    input:
        fq=OUTPUT_DIR + "/fastqs/{sample}.fastq.gz"
    output:
        counts=OUTPUT_DIR + "/readcounts/salmon/salmon_{sample}/quant.sf",
        out_dir=temp(directory(OUTPUT_DIR + "/readcounts/salmon/salmon_{sample}"))
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.salmon.quant.command_line_arguments),
    threads: CONFIG.third_party_tools.salmon.quant.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/readcounts/salmon/salmon_{sample}.log",
    shell:
        """
        salmon quant -i {params.index} {params.options} -r {input.fq} -p {threads} -o {output.out_dir}
        """

rule get_salmon_counts:
    input:
        count_dirs=expand(OUTPUT_DIR + "/readcounts/salmon/salmon_{sample}", sample=SAMPLE_NAMES),
        counts=expand(OUTPUT_DIR + "/readcounts/salmon/salmon_{sample}/quant.sf", sample=SAMPLE_NAMES)
    output:
        count_table=OUTPUT_DIR + "/readcounts/salmon/salmon_counts.csv"
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/readcounts/salmon/salmon_counts.log"
    script:
        "../scripts/get_salmon_counts.py"

ruleorder: salmon_counts_paired > salmon_counts_single
