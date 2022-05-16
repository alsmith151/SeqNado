
rule align_paired:
    input:
        fq1 = "trimmed/{sample}_R1.fastq.gz",
        fq2 = "trimmed/{sample}_R2.fastq.gz",
    params:
        index = config["aligner"]["index"],
        threads = config["pipeline"]["n_cores"],
    output:
        bam = "aligned/{sample}.bam",
    log:
        "logs/{sample}.align.log"
    shell:
        """bowtie2 -p {params.threads} -x {params.index} -1 {input.fq1} -2 {input.fq2} 2> {log} |
           samtools view -bS - > {output.bam} &&
           samtools sort -@ {params.threads} -o {output.bam}_sorted {output.bam} &&
           mv {output.bam}_sorted {output.bam} &&
           samtools index {output.bam}
        """


rule align_single:
    input:
        fq1 = "trimmed/{sample}_R0.fastq.gz",
    params:
        index = config["aligner"]["index"],
        threads = config["pipeline"]["n_cores"],
    output:
        bam = "aligned/{sample}.bam",
    log:
        "logs/{sample}.align.log"
    shell:
        """bowtie2 -p {params.threads} -x {params.index} -U {input.fq1} 2> {log} |
           samtools view -bS - > {output.bam} &&
           samtools sort -@ {params.threads} -o {output.bam}_sorted {output.bam} &&
           mv {output.bam}_sorted {output.bam} &&
           samtools index {output.bam}
        """
    




