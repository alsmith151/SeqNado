
rule align_paired:
    input:
        fq1 = "trimmed/{sample}_1.fastq.gz",
        fq2 = "trimmed/{sample}_2.fastq.gz",
    params:
        index = config["genome"]["indicies"],
    output:
        bam = "aligned/{sample}.bam",
    threads:
        config["bowtie2"]["threads"],
    log:
        "logs/align/{sample}.log"
    shell:
        """bowtie2 -p {threads} -x {params.index} -1 {input.fq1} -2 {input.fq2} 2> {log} |
           samtools view -bS - > {output.bam} &&
           samtools sort -@ {threads} -o {output.bam}_sorted {output.bam} >> {log} 2>&1 &&
           mv {output.bam}_sorted {output.bam} &&
           samtools index {output.bam}
        """

rule align_single:
    input:
        fq1 = "trimmed/{sample}.fastq.gz",
    params:
        index = config["genome"]["indicies"],
    output:
        bam = "aligned/{sample}.bam",
    threads:
        config["bowtie2"]["threads"],
    log:
        "logs/align/{sample}.log"
    shell:
        """bowtie2 -p {threads} -x {params.index} -U {input.fq1} 2> {log} |
           samtools view -bS - > {output.bam} &&
           samtools sort -@ {threads} -o {output.bam}_sorted {output.bam} &&
           mv {output.bam}_sorted {output.bam} &&
           samtools index {output.bam}
        """
    



