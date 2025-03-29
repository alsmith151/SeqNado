

rule bcftools_call_snp:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        vcf="seqnado_output/variant/{sample}.vcf.gz",
        idx="seqnado_output/variant/{sample}.vcf.gz.tbi",
    params:
        fasta=config["fasta"],
        faidx=config["fasta_index"],
    resources:
        mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt -1)}GB",
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
    threads: config["bcftools"]["threads"]
    log:
        "seqnado_output/logs/variant/{sample}.log",
    shell:"""
    bcftools mpileup --threads {threads} -Ou -f {params.fasta} {input.bam} | bcftools call --threads {threads} -mv -Oz -o {output.vcf} > {log} 2>&1
    tabix -f {input.vcf} > {output.vcf}
    """


rule bcftools_annotate:
    input:
        vcf=rules.bcftools_call_snp.output.vcf,
        idx=rules.bcftools_call_snp.output.idx,
    output:
        vcf="seqnado_output/variant/{sample}.anno.vcf.gz",
        idx="seqnado_output/variant/{sample}.anno.vcf.gz.tbi",
    params:
        dbsnp=config["snp_database"],
    threads: 16
    shell:"""
    bcftools annotate --threads 16 -c ID -a {params.dbsnp} {input.vcf} > {output.vcf}
    bcftools index -f {output.vcf}
    """