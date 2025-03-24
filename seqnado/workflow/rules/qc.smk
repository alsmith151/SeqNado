import os
from seqnado.helpers import check_options, define_time_requested, define_memory_requested

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


rule samtools_stats:
    input:
        bam="seqnado_output/aligned/raw/{sample}.bam",
    output:
        stats="seqnado_output/qc/alignment_raw/{sample}.txt",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        """samtools stats {input.bam} > {output.stats}"""


use rule samtools_stats as samtools_stats_filtered with:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
    output:
        stats="seqnado_output/qc/alignment_filtered/{sample}.txt",

def get_fastqc_files(*args, **kwargs):
    """Return a list of fastq files for a given sample name."""
    import pathlib

    fastqc_dir = pathlib.Path("seqnado_output/qc/fastqc_raw/")

    fastqc_files = []
    fq_files = pathlib.Path("seqnado_output/fastqs").glob("*.fastq.gz")
    for fq_file in fq_files:
        fastqc_file = fastqc_dir / (fq_file.stem.replace(".fastq", "") + "_fastqc.html")
        fastqc_files.append(str(fastqc_file))

    return fastqc_files


rule multiqc_raw:
    input:
        get_fastqc_files,
    output:
        "seqnado_output/qc/fastq_raw_qc.html",
    log:
        "seqnado_output/logs/multiqc_raw.log",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        "multiqc -o seqnado_output/qc seqnado_output/qc/fastqc_raw -n fastq_raw_qc.html --force > {log} 2>&1"


def get_trimmed_files(wc):
    """Return a list of fastq files for a given sample name."""
    import pathlib

    fastqc_dir = pathlib.Path("seqnado_output/qc/fastqc_trimmed/")

    fastqc_files = []
    fq_files = pathlib.Path("seqnado_output/fastqs").glob("*.fastq.gz")
    for fq_file in fq_files:
        fastqc_file = fastqc_dir / (fq_file.stem.replace(".fastq", "") + "_fastqc.html")
        fastqc_files.append(str(fastqc_file))

    return fastqc_files


rule multiqc_trimmed:
    input:
        get_trimmed_files,
    output:
        "seqnado_output/qc/fastq_trimmed_qc.html",
    log:
        "seqnado_output/logs/multiqc_trimmed.log",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        "multiqc -o seqnado_output/qc seqnado_output/qc/fastqc_trimmed -n fastq_trimmed_qc.html --force > {log} 2>&1"


rule multiqc_alignment_raw:
    input:
        expand(
            "seqnado_output/qc/alignment_raw/{sample}.txt",
            sample=SAMPLE_NAMES,
        ),
    output:
        "seqnado_output/qc/alignment_raw_qc.html",
    log:
        "seqnado_output/logs/multiqc_alignment_raw.log",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        "multiqc -o seqnado_output/qc seqnado_output/qc/alignment_raw -n alignment_raw_qc.html --force > {log} 2>&1"


rule multiqc_alignment_filtered:
    input:
        expand(
            "seqnado_output/qc/alignment_filtered/{sample}.txt",
            sample=SAMPLE_NAMES,
        ),
    output:
        "seqnado_output/qc/alignment_filtered_qc.html",
    log:
        "seqnado_output/logs/multiqc_alignment_filtered.log",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        "multiqc -o seqnado_output/qc seqnado_output/qc/alignment_filtered -n alignment_filtered_qc.html --force > {log} 2>&1"


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


ruleorder: fastqc_raw_paired > fastqc_raw_single > fastqc_trimmed_paired > fastqc_trimmed_single > samtools_stats > samtools_stats_filtered > multiqc_raw > multiqc_trimmed > multiqc_alignment_raw > multiqc_alignment_filtered > multiqc_library_complexity


rule qualimap_bamqc:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
    output:
        html="seqnado_output/qc/qualimap/bamqc_{sample}/qualimapReport.html",
    params:
        output_dir="seqnado_output/qc/qualimap/bamqc_{sample}/",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: 16
    container: "library://cchahrou/seqnado/seqnado_qc.sif:latest"
    log:"seqnado_output/logs/qualimap/bamqc_{sample}.log",
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
        html="seqnado_output/qc/qualimap/rnaseq_{sample}/qualimapReport.html",
    params:
        output_dir="seqnado_output/qc/qualimap/rnaseq_{sample}/",
        annotation=config["genome"]["gtf"],
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: 16
    container: "library://cchahrou/seqnado/seqnado_qc.sif:latest"
    log:
        "seqnado_output/logs/qualimap/rnaseq_{sample}.log",
    shell:"""
    qualimap --java-mem-size={resources.mem} rnaseq \
    -s \
    -bam {input.bam} \
    -gtf {params.annotation} \
    -outdir {params.output_dir} \
    > {log} 2>&1    
    """



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


def get_qualimap_files(wildcards):
    if ASSAY == "RNA":
        return expand(
            "seqnado_output/qc/qualimap/rnaseq_{sample}/qualimapReport.html",
            sample=SAMPLE_NAMES,
        )
    else:
        return expand(
            "seqnado_output/qc/qualimap/bamqc_{sample}/qualimapReport.html",
            sample=SAMPLE_NAMES,
        )

def get_multiqc_config():
    import importlib.resources
    import seqnado.data

    return importlib.resources.files(seqnado.data) / "multiqc_config.yaml"

rule seqnado_report:
    input:
        get_fastqc_files_all,
        expand("seqnado_output/qc/alignment_raw/{sample}.txt", sample=SAMPLE_NAMES),
        expand("seqnado_output/qc/alignment_filtered/{sample}.txt", sample=SAMPLE_NAMES),
        get_qualimap_files,
    output:
        report = "seqnado_output/seqnado_report.html",
        out_dir = temp(directory("seqnado_output/seqnado_report_data")),
    params:
        multiqc_config = get_multiqc_config(),
    log:
        "seqnado_output/logs/multiqc.log",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "library://cchahrou/seqnado/seqnado_qc.sif:latest"
    shell:"""
    multiqc -o seqnado_output seqnado_output \
    --config {params.multiqc_config} \
    --filename "seqnado_report.html" \
    --force > {log} 2>&1
    """
