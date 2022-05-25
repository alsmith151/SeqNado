rule fastqc_raw:
    input:
        fq = "fastq/{sample}.fastq.gz",
    output:
        qc = "qc/fastq_raw/{sample}_fastqc.html",
    params:
        threads = config["pipeline"]["n_cores"],
    shell:
        """
        fastqc -q -t {params.threads} --nogroup --outdir qc/fastq_raw/ {input.fq}
        """

rule multiqc_fastq_raw:
    input:
        reports = expand("qc/fastq_raw/{sample}_fastqc.html", sample=SAMPLE_NAMES_WITH_READ),
    output:
        report = "qc/fastq_qc_raw_report.html"
    shell:
        "multiqc -n fastq_qc_raw_report -o qc/ qc/fastq_raw/ --force"


rule fastqc_trimmed:
    input:
        fq = "trimmed/{sample}.fastq.gz",
    output:
        qc = "qc/fastq_trimmed/{sample}_fastqc.html",
    params:
        threads = config["pipeline"]["n_cores"],
    run:
        outdir = os.path.dirname(output.qc)
        cmd = f"""fastqc -q -t {params.threads} --nogroup --outdir {outdir} {input.fq}"""
        shell(cmd)

rule multiqc_trimmed:
    input:
        reports = expand("qc/fastq_trimmed/{sample}_fastqc.html", sample=SAMPLE_NAMES_WITH_READ),
    output:
        report = "qc/fastq_qc_trimmed_report.html"
    shell:
        "multiqc -n fastq_qc_trimmed_report.html -o qc/ qc/fastq_raw/ --force"

rule samtools_stats_raw:
    input:
        bam = "aligned/{sample}.bam"
    output:
        stats = "qc/alignment_raw/{sample}.txt"
    shell:
       """samtools stats {input.bam} > {output.stats}"""

rule multiqc_bam_raw:
    input:
        stats = expand("qc/alignment_raw/{sample_name}.txt", sample_name=SAMPLE_NAMES_NO_READ),
    output:
        report = "qc/bam_raw_qc_report.html"
    run:
        input_dir = os.path.dirname(input.stats[0])
        output_dir = os.path.dirname(output.report)
        report = os.path.basename(output.report)

        cmd = f"multiqc {input_dir} -o {output_dir} -n {report} --force"
        shell(cmd)

rule samtools_stats_filtered:
    input:
        bam = "aligned_and_filtered/{sample}.bam"
    output:
        stats = "qc/alignment_filtered/{sample}.txt"
    shell:
       """samtools stats {input.bam} > {output.stats}"""

rule multiqc_bam_filtered:
    input:
        stats = expand("qc/alignment_filtered/{sample_name}.txt", sample_name=SAMPLE_NAMES_NO_READ),
    output:
        report = "qc/bam_filtered_qc_report.html"
    run:
        input_dir = os.path.dirname(input.stats[0])
        output_dir = os.path.dirname(output.report)
        report = os.path.basename(output.report)

        cmd = f"multiqc {input_dir} -o {output_dir} -n {report} --force"
        shell(cmd)

rule multiqc_all:
    input:
        # fastq_raw_stats = expand("qc/fastq_raw/{sample}_fastqc.html", sample=SAMPLE_NAMES_WITH_READ),
        fastq_trimmed_stats = expand("qc/fastq_trimmed/{sample}_fastqc.html", sample=SAMPLE_NAMES_WITH_READ),
        # bam_raw_stats = expand("qc/alignment_raw/{sample}.txt", sample=SAMPLE_NAMES_NO_READ),
        bam_filtered_stats = expand("qc/alignment_filtered/{sample}.txt", sample=SAMPLE_NAMES_NO_READ),
    output:
        report = "qc/full_qc_report.html"
    run:
        input_dirs = [os.path.dirname(fnames[0]) for fnames in 
                     [input.fastq_trimmed_stats, input.bam_filtered_stats]
                     ]

        output_dir = os.path.dirname(output.report)
        report = os.path.basename(output.report)

        cmd = f"multiqc {' '.join(input_dirs)} -o {output_dir} -n {report} --force"
        shell(cmd)