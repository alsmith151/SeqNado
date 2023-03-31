

rule bcftools_call_snp:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
    output:
        vcf = "seqnado_output/variant/{sample}.vcf.gz",
    params:
        fasta = config["genome"]["fasta"],
        faidx = config["genome"]["fasta_index"],
    resources:
        mem_mb=1024 * 10,
    threads: config["bcftools"]["threads"]
    log:
        "seqnado_output/logs/variant/bcftools/{sample}.log",
    shell: "bcftools mpileup --threads {threads} -Ou -f {params.fasta} {input.bam} | bcftools call --threads {threads} -mv -Oz -o {output.vcf} > {log} 2>&1"

rule bcftools_split_multiallelic:
    input:
        vcf = rules.bcftools_call_snp.output.vcf,
    output:
        vcf = "seqnado_output/variant/{sample}_splitmultiallelic.vcf.gz",
    shell: "bcftools norm -m -any {input.vcf} -o {output.vcf} -Oz"

rule bcftools_filter_snp:
    input:
        vcf = rules.bcftools_split_multiallelic.output.vcf,
    output:
        vcf = "seqnado_output/variant/{sample}_filtered.vcf.gz",
    params:
        options = seqnado.utils.check_options(config["bcftools"]["filter"]),
    shell: """bcftools view -i {params.options} {input.vcf} -o {output.vcf}"""

rule index_snp:
    input:
        vcf = rules.bcftools_filter_snp.output.vcf,
    output:
        vcf = "seqnado_output/variant/{sample}_filtered.vcf.gz.tbi",
    shell: "tabix -f {input.vcf} > {output.vcf}"

rule bcftools_annotate:
    input:
        vcf = rules.bcftools_filter_snp.output.vcf,
        idx = rules.index_snp.output.vcf,
    output:
        vcf = "seqnado_output/variant/{sample}_filtered.anno.vcf.gz",
    params:
        dbsnp = config['snp_database'],
    threads: 16
    shell: "bcftools annotate --threads 16 -c ID -a {params.dbsnp} {input.vcf} > {output.vcf}"

rule bcftools_stats:
    input:
        vcf = rules.bcftools_filter_snp.output.vcf,
    output:
        stats = "seqnado_output/variant/{sample}_filtered.stats.txt",
    params:
        fasta = config["genome"]["fasta"],
    shell: "bcftools stats -F {params.fasta} -s - {input.vcf} > {output.stats}"

rule bcftools_stats_plot:
    input:
        stats = rules.bcftools_stats.output.stats,
    output:
        summary = "seqnado_output/variant/{sample}_summary.pdf",
    params:
        fasta = config["genome"]["fasta"],
        out_dir = "seqnado_output/variant/",
    shell: "plot-vcfstats -p {params.out_dir} {input.stats}"
