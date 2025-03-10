rule bwa_meth:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
    output:
        bam="seqnado_output/aligned/raw/meth/{sample}.bam"
    params:
        fasta=config["fasta"],
        options=config["bwa"]["options"]
    threads: config["bwa"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container:
        "library://asmith151/seqnado/seqnado_alignment:latest"
    log:
        "seqnado_output/logs/alignment/bwa/{sample}.log"
    shell:
        """
        bwa mem -t {threads} {params.fasta} {input.fq1} {input.fq2} | samtools view -bS - > {output.bam}
        """


rule methyldackel_extract:
    input:
        bam = "seqnado_output/aligned/raw/meth/{sample}.bam"
    output:
        bdg = "seqnado_output/methylation/{sample}_CpG.bedGraph"
    params:
        fasta=config["fasta"],
        options=config["methyldackel"]["options"]
    threads: config["methyldackel"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container:
        "library://asmith151/seqnado/seqnado_meth:latest"
    log:
        "seqnado_output/logs/methylation/methyldackel/{sample}.log"
    shell: """
        MethylDackel extract  -o {output.bdg} {params.fasta} {input.bam}
        """

