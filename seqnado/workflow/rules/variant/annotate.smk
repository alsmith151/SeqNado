
rule bcftools_annotate:
    input:
        vcf=rules.split_multiallelic.output.vcf,
        idx=rules.split_multiallelic.output.idx,
    output:
        vcf="seqnado_output/variant/{sample}.anno.vcf.gz",
        idx="seqnado_output/variant/{sample}.anno.vcf.gz.tbi",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    params:
        dbsnp=CONFIG.assay_config.snp_database,
        fasta=CONFIG.genome.fasta,
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