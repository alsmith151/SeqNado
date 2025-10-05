import os
from seqnado import Assay, AssaysWithPeakCalling, QuantificationMethod
from seqnado.helpers import define_time_requested, define_memory_requested
from seqnado.config.third_party_tools import CommandLineArguments

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
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
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
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/fastqc_raw/{sample}.log",
    shell:
        """
        fastqc -o {params.output_dir} {input} > {log} 2>&1
        """

##############################################
#                  Qualimap                  #
############################################## 

def format_qualimap_options(wildcards, options: CommandLineArguments) -> str:
    """
    Format the command line options for qualimap based on the input files and parameters.

    Mainly this removes the paired-end options if single-ended and also adds the correct options depending
    on assay type.
    """

    is_paired = INPUT_FILES.is_paired_end(wildcards.sample)
    match (is_paired, ASSAY):
        case (True, Assay.RNA):
            options = CommandLineArguments(value=options, include={"--paired", "--sorted"})
        case (True, _):
            options = CommandLineArguments(value=options, include={"--collect-overlap-pairs"})
        case (False, Assay.RNA):
            options = CommandLineArguments(value=options, exclude={"--paired", "--sorted"})
        case (False, _):
            options = CommandLineArguments(value=options, exclude={"--collect-overlap-pairs"})

    return str(options)


rule qualimap_bamqc:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
    output:
        html="seqnado_output/qc/qualimap_bamqc/{sample}/qualimapReport.html",
    params:
        output_dir="seqnado_output/qc/qualimap_bamqc/{sample}/",
        options=lambda wc: format_qualimap_options(wc, str(CONFIG.third_party_tools.qualimap.command_line_arguments)),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: 16
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
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
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
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
#            BAM filtering stats             #
##############################################

rule bam_stats:
    input: 
        sort="seqnado_output/qc/alignment_post_process/{sample}_sort.tsv",
        blacklist="seqnado_output/qc/alignment_post_process/{sample}_blacklist.tsv",
        remove_duplicates="seqnado_output/qc/alignment_post_process/{sample}_remove_duplicates.tsv",
        atac_shift="seqnado_output/qc/alignment_post_process/{sample}_atac_shift.tsv",
        filtered="seqnado_output/qc/alignment_post_process/{sample}_filter.tsv",
        final="seqnado_output/qc/alignment_post_process/{sample}_final.tsv",
    output: temp("seqnado_output/qc/alignment_post_process/{sample}_alignment_stats.tsv")
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell: """
        cat {input.sort} {input.blacklist} {input.remove_duplicates} {input.atac_shift} {input.filtered} {input.final} > {output}
    """

rule prepare_stats_report:
    input:
        expand(
            "seqnado_output/qc/alignment_post_process/{sample}_alignment_stats.tsv",
            sample=SAMPLE_NAMES,
        ),
    output:
        "seqnado_output/qc/alignment_stats.tsv",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/alignment_stats.log",
    script:
        "../scripts/alignment_stats.py"


##############################################
#               Frip Enrichment              #
##############################################

def format_frip_enrichment_options(wildcards, options: CommandLineArguments):
    is_paired = INPUT_FILES.is_paired_end(wildcards.sample)
    if is_paired:
        options.include.add("--paired")
    else:
        options.exclude.add("--paired")
    return str(options)


rule frip_enrichment:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        peak="seqnado_output/peaks/{directory}/{sample}.bed",
    output:
        pdf="seqnado_output/qc/frip_enrichment/{directory}/{sample}_frip.pdf",
        frip_count="seqnado_output/qc/frip_enrichment/{directory}/{sample}_frip.txt",
    params:
        options=lambda wc: format_frip_enrichment_options(wc, CommandLineArguments()),
    threads: 16
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
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
    single_end_assays = [name for name in SAMPLE_NAMES if not INPUT_FILES.is_paired_end(name)]
    paired_end_assays = [name for name in SAMPLE_NAMES if INPUT_FILES.is_paired_end(name)]
    fastqc_raw_paired = expand(
        "seqnado_output/qc/fastqc_raw/{sample}_{read}_fastqc.html",
        sample=paired_end_assays,
        read=[1, 2],
    )
    fastqc_raw_single = expand(
        "seqnado_output/qc/fastqc_raw/{sample}_fastqc.html",
        sample=single_end_assays,
    )
    all_qc_files = []
    for files in [fastqc_raw_paired, fastqc_raw_single]:
        if files:
            all_qc_files.extend(files)
    
    return all_qc_files


def get_fastq_screen_all(wildcards):
    single_end_assays = [name for name in SAMPLE_NAMES if not INPUT_FILES.is_paired_end(name)]
    paired_end_assays = [name for name in SAMPLE_NAMES if INPUT_FILES.is_paired_end(name)]
    fastq_screen_single = expand(
        "seqnado_output/qc/fastq_screen/{sample}_screen.txt",
        sample=single_end_assays,
    )
    fastq_screen_paired = expand(
        "seqnado_output/qc/fastq_screen/{sample}_{read}_screen.txt",
        sample=paired_end_assays,
        read=[1, 2],
    )
    all_fastq_screen_files = []
    if CONFIG.qc.run_fastq_screen:
        for files in [fastq_screen_paired, fastq_screen_single]:
            if files:
                all_fastq_screen_files.extend(files)
        
    return all_fastq_screen_files


def get_library_complexity_qc(wildcards):
    if CONFIG.qc.calculate_library_complexity:
        return expand(
            "seqnado_output/qc/library_complexity/{sample}.metrics",
            sample=SAMPLE_NAMES,
        )
    else:
        return []


def get_alignment_logs(wildcards):
    
    match ASSAY:
        case Assay.MCC:
            return []
        case Assay.RNA:
            return expand(
                "seqnado_output/aligned/star/{sample}_Log.final.out",
                sample=SAMPLE_NAMES,
            )
        case _:
            return expand(
                "seqnado_output/logs/align/{sample}.log",
                sample=SAMPLE_NAMES,
            )

def get_qualimap_files(wildcards):
    match ASSAY:
        case Assay.MCC:
            return []
        case Assay.RNA:
            return expand(
                "seqnado_output/qc/qualimap_rnaseq/{sample}/qualimapReport.html",
                sample=SAMPLE_NAMES,
            )
        case _:
            return expand(
                "seqnado_output/qc/qualimap_bamqc/{sample}/qualimapReport.html",
            sample=SAMPLE_NAMES,
        )


def get_frip_files(wildcards):
    """
    Gets the calculated FRiP (Fraction of Reads in Peaks) enrichment files.
    """
    if ASSAY in AssaysWithPeakCalling and CONFIG.assay_config.call_peaks:
        return expand(
            "seqnado_output/qc/frip_enrichment/{directory}/{sample}_frip.txt",
            sample=SAMPLE_NAMES,
            directory=[m.value for m in CONFIG.assay_config.peak_calling.method],
        )
    else:
        return []


def get_counts_files(wildcards):
    match ASSAY:
        case Assay.RNA:
            match CONFIG.assay_config.rna_quantification:
                case QuantificationMethod.FEATURECOUNTS:
                    return expand(
                        "seqnado_output/readcounts/feature_counts/read_counts.tsv",
                    )
                case QuantificationMethod.SALMON:
                    return expand(
                        "seqnado_output/readcounts/salmon/salmon_{sample}/quant.sf",
                        sample=SAMPLE_NAMES,
                    )
        case Assay.CRISPR:
            return expand(
                "seqnado_output/readcounts/feature_counts/read_counts.tsv",
            )
        case _:
            return []



def get_snp_qc(wildcards):
    if ASSAY != Assay.SNP:
        return []
    
    files = []
    
    files.extend(
        expand(
            "seqnado_output/qc/variant/{sample}.stats.txt",
            sample=SAMPLE_NAMES,
        )
    )
    if CONFIG.assay_config.annotate_snps:
        files.extend(
            expand(
                "seqnado_output/qc/variant/{sample}.anno.stats.txt",
                sample=SAMPLE_NAMES,
            )
        )
    return files

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
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:"""
    multiqc -o seqnado_output seqnado_output \
    --config {params.multiqc_config} \
    --filename "seqnado_report.html" \
    --no-data-dir \
    --force > {log} 2>&1
    """

ruleorder: fastqc_raw_paired > fastqc_raw_single > seqnado_report