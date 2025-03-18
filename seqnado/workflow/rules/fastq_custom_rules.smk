
rule deduplicate_fastq_raw:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
    output:
        fq1=temp("seqnado_output/deduplicated/{sample}/{sample}_1.fastq.gz"),
        fq2=temp("seqnado_output/deduplicated/{sample}/{sample}_2.fastq.gz"),
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/deduplication/{sample}.log",
    run:
        from mcc import mcc
        import json

        stats = mcc.deduplicate_fastq(input.fq1, output.fq1, input.fq2, output.fq2)

        with open(log[0], "w") as f:
            json.dump(stats, f)


rule flash:
    input:
        fq1="seqnado_output/deduplicated/{sample}/{sample}_1.fastq.gz",
        fq2="seqnado_output/deduplicated/{sample}/{sample}_2.fastq.gz",
    output:
        flashed=temp("seqnado_output/flashed/{sample}/{sample}.extendedFrags.fastq.gz"),
        pe1=temp("seqnado_output/flashed/{sample}/{sample}.notCombined_1.fastq.gz"),
        pe2=temp("seqnado_output/flashed/{sample}/{sample}.notCombined_2.fastq.gz"),
        hist=temp(
            "seqnado_output/flashed/{sample}/{sample}.hist"
        ),
        histogram=temp(
            "seqnado_output/flashed/{sample}/{sample}.histogram"
        ),
    params:
        outdir="seqnado_output/flashed/{sample}/{sample}",
    threads: 16
    resources:
        mem_mb=1000,
    container:  "docker://pegi3s/flash"
    log:
        "seqnado_output/logs/flash/{sample}.log",
    shell:
        """
        flash {input.fq1} {input.fq2} -o {params.outdir} -t {threads} -z --compress-prog-args pigz > {log} 2>&1
        """