import json
import shlex

def get_cutadapt_adapter_args(wildcards):
    """
    Build cutadapt adapter arguments from detected adapters.
    Falls back to config if detection fails or returns None.
    """
    adapter_file = OUTPUT_DIR + f"/resources/{wildcards.sample}_adapters.json"
    base_options = str(CONFIG.third_party_tools.cutadapt.command_line_arguments)

    try:
        with open(adapter_file, "r") as f:
            adapters = json.load(f)

        adapter_r1 = adapters.get("adapter_r1")
        adapter_r2 = adapters.get("adapter_r2")

        tokens = shlex.split(base_options)
        filtered_tokens = []
        skip_next = False
        for i, token in enumerate(tokens):
            if skip_next:
                skip_next = False
                continue
            if token in ["-g", "-G", "-a", "-A"]:
                skip_next = True
                continue
            if (
                token.startswith("-g")
                or token.startswith("-G")
                or token.startswith("-a")
                or token.startswith("-A")
            ):
                continue
            filtered_tokens.append(token)

        base_options = " ".join(filtered_tokens)

        adapter_args = ""
        if adapter_r1:
            adapter_args += f" -g '{adapter_r1}'"
        if adapter_r2:
            adapter_args += f" -G '{adapter_r2}'"

        return base_options + adapter_args

    except (FileNotFoundError, json.JSONDecodeError, KeyError):
        return base_options


rule crispr_trimming_paired:
    input:
        fq1=OUTPUT_DIR + "/fastqs/{sample}_1.fastq.gz",
        fq2=OUTPUT_DIR + "/fastqs/{sample}_2.fastq.gz",
        adapters=OUTPUT_DIR + "/resources/{sample}_adapters.json",
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
        adapters=OUTPUT_DIR + "/resources/{sample}_adapters.json",
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

ruleorder: crispr_trimming_paired > crispr_trimming_single