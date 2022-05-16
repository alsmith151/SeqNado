
rule fastqc_paired:
    input:
        f1 = "fastq/{sample}_{antibody}_R1.fastq.gz",
        f2 = "fastq/{sample}_{antibody}_R2.fastq.gz",
    output:
        r1 = "statistics/fastqc/{sample}_{antibody}_R1_fastqc.html",
        r2 = "statistics/fastqc/{sample}_{antibody}_R2_fastqc.html",
    params:
        threads = config["pipeline"]["n_cores"], 
    shell:
        """
        fastqc -q -t {params.threads} --nogroup --outdir statistics/fastqc/ {input.f1} &
        fastqc -q -t {params.threads} --nogroup --outdir statistics/fastqc/ {input.f2}
        """

rule fastqc_single:
    input:
        f1 = "fastq/{sample}_{antibody}_R0.fastq.gz",
    output:
        r1 = "statistics/fastqc/{sample}_{antibody}_R0_fastqc.html",
    params:
        threads = config["pipeline"]["n_cores"], 
    shell:
        """
        fastqc -q -t {params.threads} --nogroup --outdir statistics/fastqc/ {input.f1}
        """

rule multiqc:
    input:
        reports = expand("statistics/fastqc/{sample}_{antibody}_R{read}_fastqc.html", sample=SAMPLE_NAMES, antibody=ANTIBODIES, read=READS)
    output:
        "statistics/readqc_report.html"
    shell:
        "multiqc -f -n readqc_report -o statistics statistics/fastqc/"


ruleorder: fastqc_paired > fastqc_single
