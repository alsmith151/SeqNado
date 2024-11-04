
rule save_design:
    output:
        "seqnado_output/design.csv",
    container:
        None
    run:
        DESIGN.to_dataframe().to_csv("seqnado_output/design.csv", index=False)

rule md5sum_fastq:
    input:
        fq1="seqnado_output/fastqs/{sample}_1.fastq.gz",
        fq2="seqnado_output/fastqs/{sample}_2.fastq.gz",
    output:
        md5_1="seqnado_output/md5sums/fastq/{sample}_1.md5",
        md5_2="seqnado_output/md5sums/fastq/{sample}_2.md5",
    shell:
        """
        md5sum {input.fq1} > {output.md5_1}
        md5sum {input.fq2} > {output.md5_2}
        """

rule md5sum_fastq_combined:
    input:
        md5sums=expand("seqnado_output/md5sums/fastq/{sample}_{read}.md5", read=["1", "2"], sample=SAMPLE_NAMES),
    output:
        "seqnado_output/fastq_md5sums.txt",
    shell:
        """
        cat {input.md5sums} > {output}
        """


rule md5sum_bam:
    input:
        "seqnado_output/aligned/{sample}.bam",
    output:
        "seqnado_output/md5sums/bam/{sample}.md5",
    shell:
        """
        md5sum {input} > {output}
        """

rule md5sum_bam_combined:
    input:
        md5sums=expand("seqnado_output/md5sums/bam/{sample}.md5", sample=SAMPLE_NAMES),
    output:
        "seqnado_output/bam_md5sums.txt",
    shell:
        """
        cat {input.md5sums} > {output}
        """

rule md5sum_bigwigs:
    input:
        bigwig="seqnado_output/bigwigs/deeptools/unscaled/{sample}.bigWig"
    output:
        "seqnado_output/md5sums/bigwigs/{sample}.md5"
    shell:
        """
        md5sum {input.bigwig} > {output}
        """

rule md5sum_bigwigs_combined:
    input:
        md5sums=expand("seqnado_output/md5sums/bigwigs/{sample}.md5", sample=SAMPLE_NAMES),
    output:
        "seqnado_output/bigwig_md5sums.txt",
    shell:
        """
        cat {input.md5sums} > {output}
        """