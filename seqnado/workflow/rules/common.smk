
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
        md5_1="seqnado_output/md5sums/{sample}_1.md5",
        md5_2="seqnado_output/md5sums/{sample}_2.md5",
    shell:
        """
        md5sum {input.fq1} > {output.md5_1}
        md5sum {input.fq2} > {output.md5_2}
        """


rule md5sum_bam:
    input:
        "seqnado_output/bams/{sample}.bam",
    output:
        "seqnado_output/md5sums/{sample}.md5",
    shell:
        """
        md5sum {input} > {output}
        """