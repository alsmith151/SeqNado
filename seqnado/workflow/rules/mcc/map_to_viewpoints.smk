
rule viewpoints_to_fasta:
    input:
        bed=CONFIG.mcc_viewpoints,
        genome=CONFIG.genome.fasta,
    output:
        fasta="seqnado_output/resources/viewpoints.fa",
    log:
        "seqnado_output/logs/bed_to_fasta/viewpoints.log",
    
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output.fasta} -name 2> {log} &&
        cat {output.fasta} | sed -E 's/:+/-/g' > {output.fasta}.tmp &&
        mv {output.fasta}.tmp {output.fasta}
        """


rule fasta_index:
    input:
        fasta="seqnado_output/resources/viewpoints.fa",
    output:
        index="seqnado_output/resources/viewpoints.fa.fai",
    log:
        "seqnado_output/logs/bed_to_fasta/index.log",
    shell:
        """
        samtools faidx {input.fasta} -o {output.index}
        """

rule exclusion_regions:
    input:
        bed=CONFIG.mcc_viewpoints,
    output:
        bed="seqnado_output/resources/exclusion_regions.bed"
    log:
        "seqnado_output/logs/exclusion_regions.log"
    params:
        genome=CONFIG.genome.chromosome_sizes,
        exclusion_zone=CONFIG.assay_config.mcc.exclusion_zone
    shell:
        """
        bedtools slop -i {input.bed} -g {params.genome}  -b {params.exclusion_zone} > {output.bed}
        """


rule minimap2_to_viewpoints:
    input:
        fq="seqnado_output/flashed/{sample}/{sample}.extendedFrags.fastq.gz",
        viewpoints="seqnado_output/resources/viewpoints.fa",
    output:
        bam=temp("seqnado_output/aligned/aligned_to_viewpoints/{sample}.bam"),
    threads: 4
    resources:
        mem="4GB",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/aligned/{sample}.log",
    shell:
        """
        minimap2 -x sr -a -k 8 -w 1 --cs=long {input.viewpoints} {input.fq} 2> {log} |
        samtools view -h -b -o {output.bam} 2>> {log} &&
        samtools sort -@ 4 -o {output.bam}.sorted {output.bam} 2>> {log} &&
        mv {output.bam}.sorted {output.bam} &&
        samtools index {output.bam}
        """ 
