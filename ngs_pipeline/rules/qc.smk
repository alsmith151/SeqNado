
import os

import numpy as np

import ngs_pipeline.utils as utils

def get_reports(*args):

    reports = []

    for fq in FASTQ_FILES:

        sample_name = fq.replace(".fastq.gz", "_fastqc.html") 
        reports.append(f"qc/fastq_trimmed/{sample_name}")
       
    # Samtools reports
    for sample in SAMPLE_NAMES_NO_READ:
        reports.append(f"qc/alignment_filtered/{sample}.txt")
            
    return {"reports": reports}

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

        if not isinstance(input.reports, list):
            reports = [*input.reports]
        else:
            reports = input
        
        dirnames = [os.path.dirname(x) for x in reports]
        dirnames = list(set(dirnames))
        search_dirs = " ".join(dirnames)

        cmd = f"multiqc -n {basename} -o {outdir} {search_dirs} --force > {log} 2>&1"

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
        reports = [f"qc/fastq_raw/{sample.replace('.fastq.gz', '')}_fastqc.html" for sample in FASTQ_FILES],
    output:
        report = "qc/fastq_qc_raw_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_raw_report.log"

use rule multiqc as multiqc_fastq_trimmed with:
    input:
        reports = [f"qc/fastq_trimmed/{sample.replace('.fastq.gz', '')}_fastqc.html" for sample in FASTQ_FILES],
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
        reports = expand("qc/alignment_raw/{sample_name}.txt", sample_name=SAMPLE_NAMES_NO_READ),
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
        reports = expand("qc/alignment_filtered/{sample_name}.txt", sample_name=SAMPLE_NAMES_NO_READ),
    output:
        report = "qc/bam_filtered_qc_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_raw_report.log"

use rule multiqc as multiqc_all with:
    input:
        **get_reports(),
    output:
        report = "qc/full_qc_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_raw_report.log"