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
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
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
    log:
        "seqnado_output/logs/fastqc_raw/{sample}.log",
    shell:
        """
        fastqc -o {params.output_dir} {input} > {log} 2>&1
        """

##############################################
#                  Qualimap                  #
############################################## 

def format_qualimap_options(wildcards):
    qualimap_rnaseq_options = "--paired --sorted"
    qualimap_bamqc_options = "--collect-overlap-pairs"

    is_paired = DESIGN.query(wildcards.sample).is_paired
    if not is_paired:
        qualimap_rnaseq_options = re.sub(r"--paired", "", qualimap_rnaseq_options)
        qualimap_bamqc_options = re.sub(r"--collect-overlap-pairs", "", qualimap_bamqc_options)
    return qualimap_rnaseq_options if ASSAY == "RNA" else qualimap_bamqc_options


rule qualimap_bamqc:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
    output:
        html="seqnado_output/qc/qualimap_bamqc/{sample}/qualimapReport.html",
    params:
        output_dir="seqnado_output/qc/qualimap_bamqc/{sample}/",
        options=format_qualimap_options,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: 16
    container: "library://cchahrou/seqnado/seqnado_qc.sif:latest"
    log:"seqnado_output/logs/qualimap_bamqc/{sample}.log",
    shell:"""
    qualimap --java-mem-size={resources.mem} bamqc \
    {params.options} \
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
        options=format_qualimap_options,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: 16
    container: "library://cchahrou/seqnado/seqnado_qc.sif:latest"
    log:
        "seqnado_output/logs/qualimap_rnaseq/{sample}.log",
    shell:"""
    qualimap --java-mem-size={resources.mem} rnaseq \
    {params.options} \
    -bam {input.bam} \
    -gtf {params.annotation} \
    -outdir {params.output_dir} \
    > {log} 2>&1    
    """



##############################################
#               Frip Enrichment              #
##############################################

def format_frip_enrichment_options(wildcards):
    is_paired = DESIGN.query(wildcards.sample).is_paired
    options = "--extendReads"
    if not is_paired:
        options = re.sub(r"--extendReads", "", options)
    return options


rule frip_enrichment:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        peak="seqnado_output/peaks/{directory}/{sample}.bed",
    output:
        pdf="seqnado_output/qc/frip_enrichment/{directory}/{sample}_frip.pdf",
        frip_count="seqnado_output/qc/frip_enrichment/{directory}/{sample}_frip.txt",
    params:
        options=format_frip_enrichment_options,
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
    fastqc_raw_single = expand(
        "seqnado_output/qc/fastqc_raw/{sample}_fastqc.html",
        sample=single_end_assays,
    ),
    all_qc_files = []
    for files in [fastqc_raw_paired, fastqc_raw_single]:
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
            "seqnado_output/qc/library_complexity/{sample}.metrics",
            sample=SAMPLE_NAMES,
        )
    else:
        return []


def get_alignment_logs(wildcards):
    if ASSAY == "MCC":
        return []
    elif ASSAY == "RNA":
        return expand(
            "seqnado_output/aligned/star/{sample}_Log.final.out",
            sample=SAMPLE_NAMES,
        )
    else: 
        return expand(
            "seqnado_output/logs/align/{sample}.log",
            sample=SAMPLE_NAMES,
            )


def get_qualimap_files(wildcards):
    if ASSAY == "MCC":
        return []
    if ASSAY == "RNA":
        return expand(
            "seqnado_output/qc/qualimap_rnaseq/{sample}/qualimapReport.html",
            sample=SAMPLE_NAMES,
        )  
    else:
        return expand(
            "seqnado_output/qc/qualimap_bamqc/{sample}/qualimapReport.html",
            sample=SAMPLE_NAMES,
        )


def get_frip_files(wildcards):
    if OUTPUT.peak_calling_method:
        peak_methods = [m.value for m in OUTPUT.peak_calling_method]
        if ASSAY in ["CAT", "ATAC"] and config["call_peaks"]:
            return expand(
                "seqnado_output/qc/frip_enrichment/{directory}/{sample}_frip.txt",
                sample=SAMPLE_NAMES,
                directory=peak_methods,
            )
        if ASSAY == "ChIP" and config["call_peaks"]:
            return expand(
                "seqnado_output/qc/frip_enrichment/{directory}/{sample}_frip.txt",
                sample=SAMPLE_NAMES_IP,
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
    elif ASSAY == "CRISPR":
        return expand(
            "seqnado_output/readcounts/feature_counts/read_counts.tsv",
        )
    else:
        return []


def get_snp_qc(wildcards):
    if ASSAY == "SNP" and config["call_snps"]:
        return expand(
            "seqnado_output/qc/variant/{sample}.stats.txt",
            sample=SAMPLE_NAMES,
        )
    if ASSAY == "SNP" and config["annotate_snps"]:
        return expand(
            "seqnado_output/qc/variant/{sample}.anno.stats.txt",
        )
    else:
        return []

##############################################
#                  MultiQC                   #
##############################################

rule seqnado_report:
    input:
        get_fastqc_files_all,
        get_fastq_screen_all,
        get_alignment_logs,
        get_library_complexity_qc,
        get_qualimap_files,
        get_frip_files,
        get_counts_files,
        get_snp_qc,
    output:
        report = "seqnado_output/seqnado_report.html",
    params:
        multiqc_config = "/opt/seqnado/multiqc_config.yaml"
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

ruleorder: fastqc_raw_paired > fastqc_raw_single > seqnado_report