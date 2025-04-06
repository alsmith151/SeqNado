

rule bcftools_call_snp:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        vcf=temp("seqnado_output/variant/{sample}.raw.vcf.gz"),
        idx=temp("seqnado_output/variant/{sample}.raw.vcf.gz.tbi"),
        stats="seqnado_output/qc/variant/{sample}.stats.txt",
    params:
        fasta=config["fasta"],
    resources:
        mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt -1)}GB",
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
    threads: config["bcftools"]["threads"]
    log:
        "seqnado_output/logs/variant/{sample}.log",
    shell:"""
    bcftools mpileup --threads {threads} -Ou -f {params.fasta} {input.bam} | bcftools call --threads {threads} -mv -Oz -o {output.vcf} > {log} 2>&1 &&
    tabix -f {output.vcf} > {output.idx} &&
    bcftools stats -F {params.fasta} -s - {output.vcf} > {output.stats}
    """

rule split_multiallelic:
    input:
        vcf=rules.bcftools_call_snp.output.vcf,
        idx=rules.bcftools_call_snp.output.idx,
    output:
        vcf="seqnado_output/variant/{sample}.vcf.gz",
        idx="seqnado_output/variant/{sample}.vcf.gz.tbi",
    params:
        outdir="seqnado_output/variant/",
    resources:
        mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt -1)}GB",
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
    threads: 16
    log:
        "seqnado_output/logs/variant/{sample}_split.log",
    shell:"""
    bcftools norm --threads {threads} -m-any -Oz -o {output.vcf} {input.vcf} > {log} 2>&1 &&
    tabix -f {output.vcf} > {output.idx}
    """

rule bcftools_annotate:
    input:
        vcf=rules.split_multiallelic.output.vcf,
        idx=rules.split_multiallelic.output.idx,
    output:
        vcf="seqnado_output/variant/{sample}.anno.vcf.gz",
        idx="seqnado_output/variant/{sample}.anno.vcf.gz.tbi",
    params:
        dbsnp=config["snp_database"],
        fasta=config["fasta"],
    resources:
        mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt -1)}GB",
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
    threads: 16
    log:
        "seqnado_output/logs/variant/{sample}_anno.log",
    shell:"""
    bcftools annotate --threads {threads} -c ID -a {params.dbsnp} -Oz -o {output.vcf} {input.vcf} 2> {log} &&
    tabix -f {output.vcf} > {output.idx}
    """