from seqnado.helpers import define_time_requested, define_memory_requested

def get_trimming_input(wildcards):
    if INPUT_FILES.is_paired_end(wildcards.sample):
        return {
            "fq1": OUTPUT_DIR + f"/fastqs/{wildcards.sample}_1.fastq.gz",
            "fq2": OUTPUT_DIR + f"/fastqs/{wildcards.sample}_2.fastq.gz",
        }
    else:
        return {
            "fq": OUTPUT_DIR + f"/fastqs/{wildcards.sample}.fastq.gz",
        }

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
        fq1=OUTPUT_DIR + "/trimmed/{sample}_1.fastq.gz",
        fq2=OUTPUT_DIR + "/trimmed/{sample}_2.fastq.gz",
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
        fq=OUTPUT_DIR + "/trimmed/{sample}.fastq.gz",
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
    -a \
    {input.annotation} \
    -T \
    {threads} \
    {params.options} \
    -o \
    {output.counts} \
    {input.bam} \
    > {log} 2>&1
    """

rule mageck_count:
    input:
        fastqs=lambda wildcards: (
            expand(OUTPUT_DIR + "/trimmed/{sample}_1.fastq.gz", sample=SAMPLE_NAMES) +
            expand(OUTPUT_DIR + "/trimmed/{sample}_2.fastq.gz", sample=SAMPLE_NAMES)
            if any(INPUT_FILES.is_paired_end(s) for s in SAMPLE_NAMES)
            else expand(OUTPUT_DIR + "/trimmed/{sample}.fastq.gz", sample=SAMPLE_NAMES)
        ),
    output:
        counts=OUTPUT_DIR + "/readcounts/mageck/read_counts.txt",
    params:
        options=str(CONFIG.third_party_tools.mageck.count.command_line_arguments),
        sgrna_library=CONFIG.third_party_tools.mageck.sgrna_library if CONFIG.third_party_tools.mageck.sgrna_library else "",
        fastq_args=lambda wildcards: (
            " ".join(["-1 " + OUTPUT_DIR + f"/trimmed/{s}_1.fastq.gz -2 " + OUTPUT_DIR + f"/trimmed/{s}_2.fastq.gz" for s in SAMPLE_NAMES])
            if any(INPUT_FILES.is_paired_end(s) for s in SAMPLE_NAMES)
            else " ".join(["--fastq " + OUTPUT_DIR + f"/trimmed/{s}.fastq.gz" for s in SAMPLE_NAMES])
        ),
    threads: CONFIG.third_party_tools.mageck.count.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "quay.io/biocontainers/mageck:v0.5.9--py_0"
    log: OUTPUT_DIR + "/logs/readcounts/mageck/mageck_count.log",
    benchmark: OUTPUT_DIR + "/.benchmark/readcounts/mageck/mageck_count.tsv",
    message: "Counting guide RNAs using MAGeCK",
    shell: """
    mageck count \
    {params.fastq_args} \
    -l {params.sgrna_library} \
    -n {output.counts} \
    {params.options} \
    > {log} 2>&1
    """

rule mageck_mle:
    input:
        counts=OUTPUT_DIR + "/readcounts/mageck/read_counts.txt",
    output:
        mle_gene_summary=OUTPUT_DIR + "/readcounts/mageck/mageck_mle_gene_summary.txt",
        mle_sgrna_summary=OUTPUT_DIR + "/readcounts/mageck/mageck_mle_sgrna_summary.txt",
    params:
        options=str(CONFIG.third_party_tools.mageck.count.command_line_arguments),
        design_matrix=CONFIG.third_party_tools.mageck.design_matrix if CONFIG.third_party_tools.mageck.design_matrix else "",
    threads: CONFIG.third_party_tools.mageck.count.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "quay.io/biocontainers/mageck:v0.5.9--py_0"
    log: OUTPUT_DIR + "/logs/readcounts/mageck/mageck_mle.log",
    benchmark: OUTPUT_DIR + "/.benchmark/readcounts/mageck/mageck_mle.tsv",
    message: "Running MAGeCK MLE analysis",
    shell: """
    mageck mle \
    -k {input.counts} \
    -d {params.design_matrix} \
    -o {output.mle_gene_summary} \
    {params.options} \
    > {log} 2>&1
    """