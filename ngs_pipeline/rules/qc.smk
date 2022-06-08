rule fastqc_raw:
    input:
        fq = "fastq/{sample}.fastq.gz",
    output:
        qc = "qc/fastq_raw/{sample}_fastqc.html",
    threads:
        4,
    shell:
        """
        fastqc -q -t {threads} --nogroup --outdir qc/fastq_raw/ {input.fq}
        """

rule multiqc_fastq_raw:
    input:
        reports = expand("qc/fastq_raw/{sample}_fastqc.html", sample=SAMPLE_NAMES_WITH_READ),
    output:
        report = "qc/fastq_qc_raw_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_raw_report.log"
    shell:
        "multiqc -n fastq_qc_raw_report -o qc/ qc/fastq_raw/ --force > {log} 2>&1"


rule fastqc_trimmed:
    input:
        fq = "trimmed/{sample}.fastq.gz",
    output:
        qc = "qc/fastq_trimmed/{sample}_fastqc.html",
    threads:
        4,

    shell:
        """fastqc -q -t {threads} --nogroup --outdir qc/fastq_trimmed/ {input.fq}"""

rule multiqc_trimmed:
    input:
        reports = expand("qc/fastq_trimmed/{sample}_fastqc.html", sample=SAMPLE_NAMES_WITH_READ),
    output:
        report = "qc/fastq_qc_trimmed_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_trimmed_report.log"
    shell:
        "multiqc -n fastq_qc_trimmed_report.html -o qc/ qc/fastq_raw/ --force > {log} 2>&1"

rule samtools_stats_raw:
    input:
        bam = "aligned/{sample}.bam"
    output:
        stats = "qc/alignment_raw/{sample}.txt"
    threads:
        1,
    shell:
       """samtools stats {input.bam} > {output.stats}"""

rule multiqc_bam_raw:
    input:
        stats = expand("qc/alignment_raw/{sample_name}.txt", sample_name=SAMPLE_NAMES_NO_READ),
    output:
        report = "qc/bam_raw_qc_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_raw_report.log"
    shell:
        "multiqc qc/alignment_raw/ -o qc/ -n bam_raw_qc_report.html --force > {log} 2>&1"


rule samtools_stats_filtered:
    input:
        bam = "aligned_and_filtered/{sample}.bam"
    output:
        stats = "qc/alignment_filtered/{sample}.txt"
    threads:
        1,
    shell:
       """samtools stats {input.bam} > {output.stats}"""

rule multiqc_bam_filtered:
    input:
        stats = expand("qc/alignment_filtered/{sample_name}.txt", sample_name=SAMPLE_NAMES_NO_READ),
    output:
        report = "qc/bam_filtered_qc_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_raw_report.log"
    shell:
        "multiqc qc/alignment_filtered/ -o qc/ -n bam_filtered_qc_report.html --force > {log} 2>&1"


rule multiqc_all:
    input:
        fastq_trimmed_stats = expand("qc/fastq_trimmed/{sample}_fastqc.html", sample=SAMPLE_NAMES_WITH_READ),
        bam_filtered_stats = expand("qc/alignment_filtered/{sample}.txt", sample=SAMPLE_NAMES_NO_READ),
    output:
        report = "qc/full_qc_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_raw_report.log"
    shell:
        "multiqc {input.fastq_trimmed_stats} {input.bam_filtered_stats} -o qc/ -n full_qc_report.html --force > {log} 2>&1"