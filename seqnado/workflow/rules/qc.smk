import os
from seqnado.helpers import check_options, define_time_requested, define_memory_requested

##############################################
#                   FastQC                   #
############################################## 

rule fastqc_raw_paired:
    input:
        fq1="seqnado_output/fastqs/{sample}_1.fastq.gz",
        fq2="seqnado_output/fastqs/{sample}_2.fastq.gz",
    output:
        html1="seqnado_output/qc/fastqc_raw/{sample}_1_fastqc.html",
        html2="seqnado_output/qc/fastqc_raw/{sample}_2_fastqc.html",
        zip1="seqnado_output/qc/fastqc_raw/{sample}_1_fastqc.zip",
        zip2="seqnado_output/qc/fastqc_raw/{sample}_2_fastqc.zip",
    params:
        extra="--quiet",
        output_dir="seqnado_output/qc/fastqc_raw/",
        temp_prefix="seqnado_output/qc/fastqc_raw/{sample}",
    threads: 1
    resources:
        mem="1.5GB",
         runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/fastqc_raw/{sample}.log",
    shell:
        """
        fastqc -o {params.output_dir} {input.fq1} {input.fq2} > {log} 2>&1
        """


rule fastqc_raw_single:
    input:
        "seqnado_output/fastqs/{sample}.fastq.gz",
    output:
        html="seqnado_output/qc/fastqc_raw/{sample}_fastqc.html",
        zip="seqnado_output/qc/fastqc_raw/{sample}_fastqc.zip",
    params:
        extra="--quiet",
        output_dir="seqnado_output/qc/fastqc_raw/",
        temp_prefix="seqnado_output/qc/fastqc_raw/{sample}",
    log:
        "seqnado_output/logs/fastqc_raw/{sample}.log",
    shell:
        """
        fastqc -o {params.output_dir} {input} > {log} 2>&1
        """

use rule fastqc_raw_paired as fastqc_trimmed_paired with:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
    output:
        html1="seqnado_output/qc/fastqc_trimmed/{sample}_1_fastqc.html",
        html2="seqnado_output/qc/fastqc_trimmed/{sample}_2_fastqc.html",
        zip1="seqnado_output/qc/fastqc_trimmed/{sample}_1_fastqc.zip",
        zip2="seqnado_output/qc/fastqc_trimmed/{sample}_2_fastqc.zip",
    params:
        extra="--quiet",
        output_dir="seqnado_output/qc/fastqc_trimmed/",
        temp_prefix="seqnado_output/qc/fastqc_trimmed/{sample}",
    log:
        "seqnado_output/logs/fastqc_trimmed/{sample}.log",


use rule fastqc_raw_single as fastqc_trimmed_single with:
    input:
        "seqnado_output/trimmed/{sample}.fastq.gz",
    output:
        html="seqnado_output/qc/fastqc_trimmed/{sample}_fastqc.html",
        zip="seqnado_output/qc/fastqc_trimmed/{sample}_fastqc.zip",
    params:
        extra="--quiet",
        output_dir="seqnado_output/qc/fastqc_trimmed/",
        temp_prefix="seqnado_output/qc/fastqc_trimmed/{sample}",
    log:
        "seqnado_output/logs/fastqc_trimmed/{sample}.log",


##############################################
#            Library Complexity              #
##############################################

rule multiqc_library_complexity:
    input:
        expand(
            "seqnado_output/aligned/duplicates_removed/{sample}.metrics",
            sample=SAMPLE_NAMES,
        ),
    output:
        "seqnado_output/qc/library_complexity_qc.html",
    log:
        "seqnado_output/logs/multiqc_library_complexity.log",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        "multiqc -o seqnado_output/qc seqnado_output/aligned/duplicates_removed -n library_complexity_qc.html --force > {log} 2>&1"


##############################################
#                  Qualimap                  #
############################################## 

rule qualimap_bamqc:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
    output:
        html="seqnado_output/qc/qualimap_bamqc/{sample}/qualimapReport.html",
    params:
        output_dir="seqnado_output/qc/qualimap_bamqc/{sample}/",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: 16
    container: "library://cchahrou/seqnado/seqnado_qc.sif:latest"
    log:"seqnado_output/logs/qualimap_bamqc/{sample}.log",
    shell:"""
    qualimap --java-mem-size={resources.mem} bamqc \
    -nt {threads} \
    -bam {input.bam} \
    -outdir {params.output_dir} \
    > {log} 2>&1
    """


rule qualimap_rnaseq:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
    output:
        html="seqnado_output/qc/qualimap_rnaseq/{sample}/qualimapReport.html",
    params:
        output_dir="seqnado_output/qc/qualimap_rnaseq/{sample}/",
        annotation=config["genome"]["gtf"],
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: 16
    container: "library://cchahrou/seqnado/seqnado_qc.sif:latest"
    log:
        "seqnado_output/logs/qualimap_rnaseq/{sample}.log",
    shell:"""
    qualimap --java-mem-size={resources.mem} rnaseq \
    -s \
    -bam {input.bam} \
    -gtf {params.annotation} \
    -outdir {params.output_dir} \
    > {log} 2>&1    
    """



##############################################
#               Frip Enrichment              #
##############################################


rule frip_enrichment:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        peak="seqnado_output/peaks/{directory}/{sample}.bed",
    output:
        pdf="seqnado_output/qc/frip_enrichment/{directory}/{sample}_frip.pdf",
        frip_count="seqnado_output/qc/frip_enrichment/{directory}/{sample}_frip.txt",
    threads: 16
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:"seqnado_output/logs/frip_enrichment/{directory}/{sample}.log",
    shell:
        """
        plotEnrichment -p {threads} \
        --bamfiles {input.bam} \
        --BED {input.peak} \
        --outRawCounts {output.frip_count} \
        --plotFile {output.pdf} \
        > {log} 2>&1
        """


##############################################
#                Gather Stats                #
##############################################

def get_fastqc_files_all(wildcards):
    single_end_assays = [name for name in SAMPLE_NAMES if DESIGN.query(name).is_paired == False]
    paired_end_assays = [name for name in SAMPLE_NAMES if DESIGN.query(name).is_paired == True]
    fastqc_raw_paired = expand(
        "seqnado_output/qc/fastqc_raw/{sample}_{read}_fastqc.html",
        sample=paired_end_assays,
        read=[1, 2],
    ),
    fastqc_trimmed_paired = expand(
        "seqnado_output/qc/fastqc_trimmed/{sample}_{read}_fastqc.html",
        sample=paired_end_assays,
        read=[1, 2],
    ),
    fastqc_raw_single = expand(
        "seqnado_output/qc/fastqc_raw/{sample}_fastqc.html",
        sample=single_end_assays,
    ),
    fastqc_trimmed_single = expand(
        "seqnado_output/qc/fastqc_trimmed/{sample}_fastqc.html",
        sample=single_end_assays,
    ),
    
    all_qc_files = []
    for files in [fastqc_raw_paired, fastqc_trimmed_paired, fastqc_raw_single, fastqc_trimmed_single]:
        if files:
            all_qc_files.extend(*files)
    
    return all_qc_files

def get_fastq_screen_all(wildcards):
    single_end_assays = [name for name in SAMPLE_NAMES if DESIGN.query(name).is_paired == False]
    paired_end_assays = [name for name in SAMPLE_NAMES if DESIGN.query(name).is_paired == True]
    fastq_screen_single = expand(
        "seqnado_output/qc/fastq_screen/{sample}_screen.txt",
        sample=single_end_assays,
    ),
    fastq_screen_paired = expand(
        "seqnado_output/qc/fastq_screen/{sample}_{read}_screen.txt",
        sample=paired_end_assays,
        read=[1, 2],
    ),
    all_fastq_screen_files = []
    if config["fastq_screen"]:
        for files in [fastq_screen_paired, fastq_screen_single]:
            if files:
                all_fastq_screen_files.extend(*files)
        
    return all_fastq_screen_files


def get_library_complexity_qc(wildcards):
    if config["library_complexity"]:
        return expand(
            "seqnado_output/aligned/duplicates_removed/{sample}.metrics",
            sample=SAMPLE_NAMES,
        )

def get_alignment_logs(wildcards):
    return expand(
        "seqnado_output/aligned/star/{sample}_Log.final.out",
        sample=SAMPLE_NAMES,
    ) if ASSAY == "RNA" else expand(
        "seqnado_output/logs/align/{sample}.log",
        sample=SAMPLE_NAMES,
    )


rule prepare_stats_report:
    input:
        expand(
            "seqnado_output/qc/alignment_post_process/{sample}_alignment_stats.tsv",
            sample=SAMPLE_NAMES,
        ),
    output:
        "seqnado_output/qc/alignment_post_process/alignment_stats.tsv",
    log:
        "seqnado_output/logs/alignment_stats.log",
    script:
        "../scripts/alignment_stats.py"


def get_qualimap_files(wildcards):
    return expand(
        "seqnado_output/qc/qualimap_rnaseq/{sample}/qualimapReport.html",
        sample=SAMPLE_NAMES,
    ) if ASSAY == "RNA" else expand(
        "seqnado_output/qc/qualimap_bamqc/{sample}/qualimapReport.html",
        sample=SAMPLE_NAMES,
    )


def get_frip_files(wildcards):
    import glob
    if ASSAY != "RNA" and config["call_peaks"]:
        peak_methods = OUTPUT.peak_calling_method
        return expand(
            "seqnado_output/qc/frip_enrichment/{directory}/{sample}_frip.txt",
            sample=SAMPLE_NAMES,
            directory=peak_methods,
        )
    else:
        return []


def get_counts_files(wildcards):
    if ASSAY == "RNA" and config["rna_quantification"] == "featurecounts":
        return expand(
            "seqnado_output/readcounts/feature_counts/read_counts.tsv",
        ) 
    elif ASSAY == "RNA" and config["rna_quantification"] == "salmon":
        return expand(
            "seqnado_output/readcounts/salmon/salmon_{sample}/quant.sf",
            sample=SAMPLE_NAMES,
        ) 
    else:
        return []


##############################################
#                  MultiQC                   #
##############################################

def get_multiqc_config():
    import importlib.resources
    import seqnado.data
    import pathlib

    return pathlib.Path(importlib.resources.files(seqnado.data) / "multiqc_config.yaml").absolute().resolve()


rule seqnado_report:
    input:
        get_fastqc_files_all,
        get_fastq_screen_all,
        get_alignment_logs,
        get_library_complexity_qc,
        get_qualimap_files,
        "seqnado_output/qc/alignment_post_process/alignment_stats.tsv",
        get_frip_files,
        get_counts_files,
    output:
        report = "seqnado_output/seqnado_report.html",
    params: multiqc_config = get_multiqc_config(),
    log: "seqnado_output/logs/seqnado_report.log",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "library://cchahrou/seqnado/seqnado_qc.sif:latest"
    shell:"""
    multiqc -o seqnado_output seqnado_output \
    --config {params.multiqc_config} \
    --filename "seqnado_report.html" \
    --no-data-dir \
    --force > {log} 2>&1
    """

ruleorder: fastqc_raw_paired > fastqc_raw_single > fastqc_trimmed_paired > fastqc_trimmed_single > multiqc_library_complexity > seqnado_report