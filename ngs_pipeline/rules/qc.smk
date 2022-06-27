
import os

import numpy as np

import ngs_pipeline.utils as utils

def get_reports(*args):

    reports = []

    for sample in SAMPLE_NAMES_NO_READ:
        # Fastqc reports
        if df_samples_paired.query(f"ip == '{sample}'")["paired_or_single"].values[0] == "paired":
           reports.append(f"qc/fastq_trimmed/{sample}_1_fastqc.html")
           reports.append(f"qc/fastq_trimmed/{sample}_2_fastqc.html")
        else:
            reports.append(f"qc/fastq_trimmed/{sample}_fastqc.html")

        # Samtools reports
        reports.append(f"qc/alignment_filtered/{sample}.txt")
    
    return reports


# rule multiqc_fastq_raw:
#     input:
#         reports = expand("qc/fastq_raw/{sample}_fastqc.html", sample=SAMPLE_NAMES_WITH_READ),
#     output:
#         report = "qc/fastq_qc_raw_report.html"
#     threads:
#         1,
#     log:
#         "logs/qc/fastq_qc_raw_report.log"
#     shell:
#         "multiqc -n fastq_qc_raw_report -o qc/ qc/fastq_raw/ --force > {log} 2>&1"


# rule fastqc_trimmed:
#     input:
#         fq = "trimmed/{sample}.fastq.gz",
#     output:
#         qc = "qc/fastq_trimmed/{sample}_fastqc.html",
#     threads:
#         4,

#     shell:
#         """fastqc -q -t {threads} --nogroup --outdir qc/fastq_trimmed/ {input.fq}"""

# rule multiqc_trimmed:
#     input:
#         reports = expand("qc/fastq_trimmed/{sample}_fastqc.html", sample=SAMPLE_NAMES_WITH_READ),
#     output:
#         report = "qc/fastq_qc_trimmed_report.html"
#     threads:
#         1,
#     log:
#         "logs/qc/fastq_qc_trimmed_report.log"
#     shell:
#         "multiqc -n fastq_qc_trimmed_report.html -o qc/ qc/fastq_raw/ --force > {log} 2>&1"

# rule samtools_stats_raw:
#     input:
#         bam = "aligned/{sample}.bam"
#     output:
#         stats = "qc/alignment_raw/{sample}.txt"
#     threads:
#         1,
#     shell:
#        """samtools stats {input.bam} > {output.stats}"""

# rule multiqc_bam_raw:
#     input:
#         stats = expand("qc/alignment_raw/{sample_name}.txt", sample_name=SAMPLE_NAMES_NO_READ),
#     output:
#         report = "qc/bam_raw_qc_report.html"
#     threads:
#         1,
#     log:
#         "logs/qc/fastq_qc_raw_report.log"
#     shell:
#         "multiqc qc/alignment_raw/ -o qc/ -n bam_raw_qc_report.html --force > {log} 2>&1"


# rule samtools_stats_filtered:
#     input:
#         bam = "aligned_and_filtered/{sample}.bam"
#     output:
#         stats = "qc/alignment_filtered/{sample}.txt"
#     threads:
#         1,
#     shell:
#        """samtools stats {input.bam} > {output.stats}"""

# rule multiqc_bam_filtered:
#     input:
#         stats = expand("qc/alignment_filtered/{sample_name}.txt", sample_name=SAMPLE_NAMES_NO_READ),
#     output:
#         report = "qc/bam_filtered_qc_report.html"
#     threads:
#         1,
#     log:
#         "logs/qc/fastq_qc_raw_report.log"
#     shell:
#         "multiqc qc/alignment_filtered/ -o qc/ -n bam_filtered_qc_report.html --force > {log} 2>&1"


# rule multiqc_all:
#     input:
#         fastq_trimmed_stats = expand("qc/fastq_trimmed/{sample}_fastqc.html", sample=SAMPLE_NAMES_WITH_READ),
#         bam_filtered_stats = expand("qc/alignment_filtered/{sample}.txt", sample=SAMPLE_NAMES_NO_READ),
#     output:
#         report = "qc/full_qc_report.html"
#     threads:
#         1,
#     log:
#         "logs/qc/fastq_qc_raw_report.log"
#     shell:
#         "multiqc {input.fastq_trimmed_stats} {input.bam_filtered_stats} -o qc/ -n full_qc_report.html --force > {log} 2>&1"


rule fastqc:
    input:
        fq = "FASTQ",
    output:
        qc = "OUTPUT",
    params:
        outdir = "OUTDIR",
    threads:
        4,
    shell:
        """
        fastqc -q -t {threads} --nogroup --outdir {params.outdir} {input.fq}
        """

rule multiqc:
    input:
        reports = "INPUTS",
    output:
        report = "OUTPUT",
    threads:
        1,
    log:
        "LOG"
    run:

        outdir = os.path.dirname(output.report)
        basename = os.path.basename(output.report)

        dirnames = [os.path.dirname(x) for x in input.reports]
        dirnames = list(set(dirnames))
        search_dirs = " ".join(dirnames)

        cmd = f"multiqc -n {basename} -o {outdir} {search_dirs} --force > {log}"

        if workflow.use_singularity:
            cmd = utils.get_singularity_command(command=cmd,
                                                workflow=workflow,)
        shell(cmd)


use rule fastqc as fastqc_raw with:
    input:
        fq = "fastq/{sample}.fastq.gz",
    output:
        qc = "qc/fastq_raw/{sample}_fastqc.html",
    params:
        outdir = "qc/fastq_raw"

use rule fastqc as fastqc_trimmed with:
    input:
        fq = "trimmed/{sample}.fastq.gz",
    output:
        qc = "qc/fastq_trimmed/{sample}_fastqc.html",
    params:
        outdir = "qc/fastq_trimmed"

use rule multiqc as multiqc_fastq_raw with:
    input:
        reports = expand("qc/fastq_raw/{sample}_fastqc.html", sample=SAMPLE_NAMES_WITH_READ),
    output:
        report = "qc/fastq_qc_raw_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_raw_report.log"

use rule multiqc as multiqc_fastq_trimmed with:
    input:
        reports = expand("qc/fastq_trimmed/{sample}_fastqc.html", sample=SAMPLE_NAMES_WITH_READ),
    output:
        report = "qc/fastq_qc_trimmed_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_trimmed_report.log"

rule samtools_stats:
    input:
        bam = "aligned/{sample}.bam"
    output:
        stats = "qc/alignment_raw/{sample}.txt"
    threads:
        1,
    shell:
       """samtools stats {input.bam} > {output.stats}"""

use rule multiqc as multiqc_bam_raw with:
    input:
        stats = expand("qc/alignment_raw/{sample_name}.txt", sample_name=SAMPLE_NAMES_NO_READ),
    output:
        report = "qc/bam_raw_qc_report.html"
    threads:
        1,
    log:
        "logs/qc/bam_qc_raw_report.log"

use rule samtools_stats as samtools_stats_filtered with:
    input:
        bam = "aligned_and_filtered/{sample}.bam"
    output:
        stats = "qc/alignment_filtered/{sample}.txt"


use rule multiqc as multiqc_bam_filtered with:
    input:
        stats = expand("qc/alignment_filtered/{sample_name}.txt", sample_name=SAMPLE_NAMES_NO_READ),
    output:
        report = "qc/bam_filtered_qc_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_raw_report.log"

use rule multiqc as multiqc_all with:
    input:
        reports = get_reports,
    output:
        report = "qc/full_qc_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_raw_report.log"