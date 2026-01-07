from seqnado.helpers import (
    define_time_requested,
    define_memory_requested,
    should_trim_fastqs,
    get_alignment_input_paired,
    get_alignment_input_single,
)
import json

# Check if trimming is enabled
TRIMMING_ENABLED = should_trim_fastqs(CONFIG)

# Check cutadapt mode to determine if we should auto-detect adapters
CUTADAPT_MODE = getattr(CONFIG.third_party_tools.cutadapt, 'mode', None)
AUTO_DETECT_ADAPTERS = CUTADAPT_MODE == 'crispr'

if TRIMMING_ENABLED and AUTO_DETECT_ADAPTERS:
    # Rule to detect adapters from a representative sample
    rule detect_crispr_adapters_paired:
        input:
            fq1=OUTPUT_DIR + "/fastqs/{sample}_1.fastq.gz",
            fq2=OUTPUT_DIR + "/fastqs/{sample}_2.fastq.gz",
        output:
            adapters=OUTPUT_DIR + "/qc/adapters/{sample}_adapters.json",
        params:
            n_reads=10000,
            min_length=6,
            max_length=100,
            min_frequency=0.7,
        threads: 1
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/adapter_detection/{sample}.log",
        benchmark: OUTPUT_DIR + "/.benchmark/adapter_detection/{sample}.tsv",
        message: "Detecting adapters for CRISPR sample {wildcards.sample} (paired-end)",
        script:
            "../scripts/detect_crispr_adapters.py"

    rule detect_crispr_adapters_single:
        input:
            fq1=OUTPUT_DIR + "/fastqs/{sample}.fastq.gz",
        output:
            adapters=OUTPUT_DIR + "/qc/adapters/{sample}_adapters.json",
        params:
            n_reads=10000,
            min_length=6,
            max_length=100,
            min_frequency=0.7,
        threads: 1
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/adapter_detection/{sample}.log",
        benchmark: OUTPUT_DIR + "/.benchmark/adapter_detection/{sample}.tsv",
        message: "Detecting adapters for CRISPR sample {wildcards.sample} (single-end)",
        script:
            "../scripts/detect_crispr_adapters.py"

    ruleorder: detect_crispr_adapters_paired > detect_crispr_adapters_single


    def get_cutadapt_adapter_args(wildcards):
        """
        Build cutadapt adapter arguments from detected adapters.
        Falls back to config if detection fails or returns None.
        """
        adapter_file = OUTPUT_DIR + f"/qc/adapters/{wildcards.sample}_adapters.json"

        # Get base options from config (without adapter specifications)
        base_options = str(CONFIG.third_party_tools.cutadapt.command_line_arguments)

        # Try to load detected adapters
        try:
            with open(adapter_file, 'r') as f:
                adapters = json.load(f)

            adapter_r1 = adapters.get('adapter_r1')
            adapter_r2 = adapters.get('adapter_r2')

            # Remove any existing -g/-G/-a/-A flags from base_options
            import shlex
            tokens = shlex.split(base_options)
            filtered_tokens = []
            skip_next = False
            for i, token in enumerate(tokens):
                if skip_next:
                    skip_next = False
                    continue
                if token in ['-g', '-G', '-a', '-A']:
                    skip_next = True
                    continue
                if token.startswith('-g') or token.startswith('-G') or token.startswith('-a') or token.startswith('-A'):
                    continue
                filtered_tokens.append(token)

            base_options = ' '.join(filtered_tokens)

            # Add detected adapters
            adapter_args = ""
            if adapter_r1:
                adapter_args += f" -g '{adapter_r1}'"
            if adapter_r2:
                adapter_args += f" -G '{adapter_r2}'"

            return base_options + adapter_args

        except (FileNotFoundError, json.JSONDecodeError, KeyError):
            # Fall back to config if detection failed
            return base_options


    rule crispr_trimming_paired:
        input:
            fq1=OUTPUT_DIR + "/fastqs/{sample}_1.fastq.gz",
            fq2=OUTPUT_DIR + "/fastqs/{sample}_2.fastq.gz",
            adapters=OUTPUT_DIR + "/qc/adapters/{sample}_adapters.json",
        output:
            trimmed1=temp(OUTPUT_DIR + "/trimmed/{sample}_1.fastq.gz"),
            trimmed2=temp(OUTPUT_DIR + "/trimmed/{sample}_2.fastq.gz"),
        params:
            options=get_cutadapt_adapter_args,
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
            adapters=OUTPUT_DIR + "/qc/adapters/{sample}_adapters.json",
        output:
            trimmed=temp(OUTPUT_DIR + "/trimmed/{sample}.fastq.gz"),
        params:
            options=get_cutadapt_adapter_args,
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

elif TRIMMING_ENABLED:
    # Use static adapters from config (original behavior)
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
