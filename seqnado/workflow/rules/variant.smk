if config.get("call_snps"):

    rule bcftools_call_snp:
        input:
            bam="seqnado_output/aligned/{sample}.bam",
            bai="seqnado_output/aligned/{sample}.bam.bai",
        output:
            vcf="seqnado_output/variant/{sample}.vcf.gz",
        params:
            fasta=config["fasta"],
            faidx=config["fasta_index"],
        resources:
            mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt -1)}GB",
            runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
        threads: config["bcftools"]["threads"]
        log:
            "seqnado_output/logs/variant/{sample}.log",
        shell:
            "bcftools mpileup --threads {threads} -Ou -f {params.fasta} {input.bam} | bcftools call --threads {threads} -mv -Oz -o {output.vcf} > {log} 2>&1"

    rule index_snp:
        input:
            vcf=rules.bcftools_call_snp.output.vcf,
        output:
            vcf="seqnado_output/variant/{sample}.vcf.gz.tbi",
        shell:
            "tabix -f {input.vcf} > {output.vcf}"


    rule bcftools_annotate:
        input:
            vcf=rules.bcftools_call_snp.output.vcf,
            idx=rules.index_snp.output.vcf,
        output:
            vcf="seqnado_output/variant/{sample}.anno.vcf.gz",
        params:
            dbsnp=config["snp_database"],
        threads: 16
        shell:
            "bcftools annotate --threads 16 -c ID -a {params.dbsnp} {input.vcf} > {output.vcf}"
