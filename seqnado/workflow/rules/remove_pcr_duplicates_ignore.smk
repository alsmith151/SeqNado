rule ignore_duplicates:
    input:
        bam="aligned/{sample}.bam",
        index="aligned/{sample}.bam.bai",
    output:
        bam="aligned_and_filtered/{sample}.bam",
        log="logs/duplicate_removal/not_removed/{sample}.log",
    log:
        "logs/duplicate_removal/not_removed/{sample}.log",
    shell:
        """
        ln -s $(realpath {input.bam}) {output.bam} &&
        ln -s {input.bam}.bai {output.bam}.bai
        """

localrules: ignore_duplicates
