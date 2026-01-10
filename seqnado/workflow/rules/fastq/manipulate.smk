rule deduplicate_fastq_raw:
    input:
        fq1=OUTPUT_DIR + "/trimmed/{sample}_1.fastq.gz",
        fq2=OUTPUT_DIR + "/trimmed/{sample}_2.fastq.gz",
    output:
        fq1=temp(OUTPUT_DIR + "/deduplicated/{sample}/{sample}_1.fastq.gz"),
        fq2=temp(OUTPUT_DIR + "/deduplicated/{sample}/{sample}_2.fastq.gz"),
    threads: 1
    resources:
        mem="1GB",
    container:  "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/deduplication/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/deduplication/{sample}.tsv",
    message: "Deduplicating reads for sample {wildcards.sample}",
    script:
        "../../scripts/deduplicate_fastq.py"


rule flash:
    input:
        fq1=OUTPUT_DIR + "/deduplicated/{sample}/{sample}_1.fastq.gz",
        fq2=OUTPUT_DIR + "/deduplicated/{sample}/{sample}_2.fastq.gz",
    output:
        flashed=temp(OUTPUT_DIR + "/flashed/{sample}/{sample}.extendedFrags.fastq.gz"),
        pe1=temp(OUTPUT_DIR + "/flashed/{sample}/{sample}.notCombined_1.fastq.gz"),
        pe2=temp(OUTPUT_DIR + "/flashed/{sample}/{sample}.notCombined_2.fastq.gz"),
        hist=temp(
            OUTPUT_DIR + "/flashed/{sample}/{sample}.hist"
        ),
        histogram=temp(
            OUTPUT_DIR + "/flashed/{sample}/{sample}.histogram"
        ),
    params:
        outdir=OUTPUT_DIR + "/flashed/{sample}/{sample}",
    threads: 16
    resources:
        mem_mb=1000,
    container:  "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/flash/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/flash/{sample}.tsv",
    message: "Merging overlapping paired-end reads for sample {wildcards.sample} using FLASH",
    shell:
        """
        flash {input.fq1} {input.fq2} -o {params.outdir} -t {threads} -z --compress-prog-args pigz > {log} 2>&1
        """
