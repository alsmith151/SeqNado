import seqnado.utils as utils


rule align_paired:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
    params:
        index=config["genome"]["indices"],
        options=utils.check_options(config["bowtie2"]["options"]),
    output:
        bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
    threads: config["bowtie2"]["threads"]
    resources:
        time=lambda wildcards, attempt: "0-{hours}:00:00".format(
            hours=4 * 2 ** (attempt - 1)
        ),
        mem_mb=4000,
    log:
        "seqnado_output/logs/align/{sample}.log",
    shell:
        """bowtie2 -p {threads} -x {params.index} -1 {input.fq1} -2 {input.fq2} {params.options} 2> {log} |
           samtools view -bS - > {output.bam} &&
           samtools sort -@ {threads} -o {output.bam}_sorted {output.bam} >> {log} 2>&1 &&
           mv {output.bam}_sorted {output.bam}
        """


rule align_single:
    input:
        fq1="seqnado_output/trimmed/{sample}.fastq.gz",
    params:
        index=config["genome"]["indices"],
        options=utils.check_options(config["bowtie2"]["options"]),
    output:
        bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
    resources:
        time=lambda wildcards, attempt: "0-{hours}:00:00".format(
            hours=4 * 2 ** (attempt - 1)
        ),
        mem_mb=4000,
    threads: config["bowtie2"]["threads"]
    log:
        "seqnado_output/logs/align/{sample}.log",
    shell:
        """bowtie2 -p {threads} -x {params.index} -U {input.fq1} {params.options} 2> {log} |
            samtools view -bS - > {output.bam} &&
            samtools sort -@ {threads} -o {output.bam}_sorted {output.bam} &&
            mv {output.bam}_sorted {output.bam}
        """


ruleorder: align_paired > align_single
