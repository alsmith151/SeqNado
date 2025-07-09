from seqnado.helpers import check_options, define_time_requested, define_memory_requested



rule align_paired:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
    params:
        index=config["genome"]["index"],
        options=check_options(config["bowtie2"]["options"]),
        rg="--rg-id {sample} --rg SM:{sample}",
    output:
        bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
    threads: config["bowtie2"]["threads"]
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/align/{sample}.log",
    shell:
        """bowtie2 -p {threads} -x {params.index} -1 {input.fq1} -2 {input.fq2} {params.rg} {params.options} 2> {log} |
           samtools view -bS - > {output.bam}
        """


use rule align_paired as align_single with:
    input:
        fq1="seqnado_output/trimmed/{sample}.fastq.gz",
    output:
        bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
    shell:
        """bowtie2 -p {threads} -x {params.index} -U {input.fq1} {params.rg} {params.options} 2> {log} |
            samtools view -bS - > {output.bam}
        """


ruleorder: align_paired > align_single
