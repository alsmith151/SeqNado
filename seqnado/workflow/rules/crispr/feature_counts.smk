from seqnado.helpers import (
    define_time_requested, 
    define_memory_requested,
    should_trim_fastqs,
    get_alignment_input_paired,
    get_alignment_input_single,
    get_trimmed_fastq_paired,
    get_trimmed_fastq_single,
    get_raw_fastq_paired,
    get_raw_fastq_single,
)

# Check if trimming is enabled
TRIMMING_ENABLED = should_trim_fastqs(CONFIG)

if TRIMMING_ENABLED:
    rule crispr_trimming_paired:
        input:
            fq1=OUTPUT_DIR + "/fastqs/{sample}_1.fastq.gz",
            fq2=OUTPUT_DIR + "/fastqs/{sample}_2.fastq.gz",
        output:
            trimmed1=temp(OUTPUT_DIR + "/trimmed/{sample}_1.fastq.gz"),
            trimmed2=temp(OUTPUT_DIR + "/trimmed/{sample}_2.fastq.gz"),
        params:
            options=str(CONFIG.third_party_tools.cutadapt.command_line_arguments),
            trim_dir=OUTPUT_DIR + "/trimmed",
        threads: CONFIG.third_party_tools.cutadapt.threads
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/trimming/{sample}.log",
        benchmark: OUTPUT_DIR + "/.benchmark/trimming/{sample}.tsv",
        message: "Trimming adapters from CRISPR sample {wildcards.sample} using cutadapt (paired-end)",
        shell: """
        cutadapt {params.options} -o {output.trimmed1} -p {output.trimmed2} {input.fq1} {input.fq2} > {log} 2>&1
        """

    rule crispr_trimming_single:
        input:
            fq=OUTPUT_DIR + "/fastqs/{sample}.fastq.gz",
        output:
            trimmed=temp(OUTPUT_DIR + "/trimmed/{sample}.fastq.gz"),
        params:
            options=str(CONFIG.third_party_tools.cutadapt.command_line_arguments),
            trim_dir=OUTPUT_DIR + "/trimmed",
        threads: CONFIG.third_party_tools.cutadapt.threads
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/trimming/{sample}.log",
        benchmark: OUTPUT_DIR + "/.benchmark/trimming/{sample}.tsv",
        message: "Trimming adapters from CRISPR sample {wildcards.sample} using cutadapt (single-end)",
        shell: """
        cutadapt {params.options} -o {output.trimmed} {input.fq} > {log} 2>&1
        """

rule align_crispr_paired:
    input:
        fq1=lambda wildcards: get_alignment_input_paired(wildcards, OUTPUT_DIR, TRIMMING_ENABLED)["fq1"],
        fq2=lambda wildcards: get_alignment_input_paired(wildcards, OUTPUT_DIR, TRIMMING_ENABLED)["fq2"],
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.align.command_line_arguments),
    output:
        bam_unsorted=temp(OUTPUT_DIR + "/aligned/{sample}.unsorted.bam"),
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: CONFIG.third_party_tools.bowtie2.align.threads
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/align/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}.tsv",
    message: "Aligning CRISPR sample {wildcards.sample} using Bowtie2 (paired-end)",
    shell: """
    bowtie2 -p {threads} -x {params.index} -1 {input.fq1} -2 {input.fq2} {params.options} 2> {log} |
    samtools view -bS - > {output.bam_unsorted} &&
    samtools sort -@ {threads} {output.bam_unsorted} -o {output.bam}
    """

rule align_crispr_single:
    input:
        fq=lambda wildcards: get_alignment_input_single(wildcards, OUTPUT_DIR, TRIMMING_ENABLED),
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.align.command_line_arguments),
    output:
        bam_unsorted=temp(OUTPUT_DIR + "/aligned/{sample}.unsorted.bam"),
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: CONFIG.third_party_tools.bowtie2.align.threads
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/align/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}.tsv",
    message: "Aligning CRISPR sample {wildcards.sample} using Bowtie2 (single-end)",
    shell: """
    bowtie2 -p {threads} -x {params.index} -U {input.fq} {params.options} 2> {log} |
    samtools view -bS - > {output.bam_unsorted} &&
    samtools sort -@ {threads} {output.bam_unsorted} -o {output.bam}
    """

ruleorder: align_crispr_paired > align_crispr_single

if TRIMMING_ENABLED:
    ruleorder: crispr_trimming_paired > crispr_trimming_single

rule index_crispr_bam:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
    output:
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
    threads: CONFIG.third_party_tools.samtools.index.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/align/{sample}_index.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}_index.tsv",
    message: "Indexing BAM file for CRISPR sample {wildcards.sample}",
    shell: """
    samtools index {input.bam} {output.bai} > {log} 2>&1
    """

rule feature_counts:
    input:
        bam=expand(OUTPUT_DIR + "/aligned/{sample}.bam", sample=SAMPLE_NAMES),
        bai=expand(OUTPUT_DIR + "/aligned/{sample}.bam.bai", sample=SAMPLE_NAMES),
        annotation=config["genome"]["gtf"],
    output:
        counts=OUTPUT_DIR + "/readcounts/feature_counts/read_counts.tsv",
    params:
        options=str(CONFIG.third_party_tools.subread.feature_counts.command_line_arguments),
        paired= "-p --countReadPairs" if INPUT_FILES.is_paired_end(SAMPLE_NAMES[0]) else "",
        annotation_format= "SAF" if config["genome"]["gtf"].endswith(".saf") else "GTF",
    threads: CONFIG.third_party_tools.subread.feature_counts.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/readcounts/featurecounts/featurecounts.log",   
    benchmark: OUTPUT_DIR + "/.benchmark/readcounts/featurecounts/featurecounts.tsv",
    message: "Counting features using featureCounts",
    shell: """
    featureCounts \
    -a {input.annotation} \
    -F {params.annotation_format} \
    -T {threads} \
    {params.options} \
    {params.paired} \
    --donotsort \
    -o {output.counts} \
    {input.bam} \
    > {log} 2>&1
    """
