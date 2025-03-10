from seqnado.helpers import check_options, define_time_requested, define_memory_requested


rule trimgalore_paired:
    # Trim reads using trimgalore
    input:
        fq1="seqnado_output/fastqs/{sample}_1.fastq.gz",
        fq2="seqnado_output/fastqs/{sample}_2.fastq.gz",
    output:
        trimmed1=temp("seqnado_output/trimmed/{sample}_1.fastq.gz"),
        trimmed2=temp("seqnado_output/trimmed/{sample}_2.fastq.gz"),
    threads: config["trim_galore"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    params:
        options=check_options(config["trim_galore"]["options"]),
        trim_dir="seqnado_output/trimmed",
    log:
        "seqnado_output/logs/trimming/{sample}.log",
    shell:
        """
        trim_galore --cores {threads} {params.options} --basename {wildcards.sample} --paired --output_dir {params.trim_dir} {input.fq1} {input.fq2} >> {log} 2>&1 &&
        mv {params.trim_dir}/{wildcards.sample}_val_1.fq.gz {output.trimmed1} &&
        mv {params.trim_dir}/{wildcards.sample}_val_2.fq.gz {output.trimmed2}
        """


rule trimgalore_single:
    # Trim reads using trimgalore
    input:
        fq="seqnado_output/fastqs/{sample}.fastq.gz",
    output:
        trimmed=temp("seqnado_output/trimmed/{sample}.fastq.gz"),
    threads: config["trim_galore"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    params:
        options=check_options(config["trim_galore"]["options"]),
        trim_dir="seqnado_output/trimmed",
    log:
        "seqnado_output/logs/trimming/{sample}.log",
    shell:
        """
        trim_galore --cores {threads} {params.options} --basename {wildcards.sample} --output_dir {params.trim_dir} {input.fq} >> {log} 2>&1 &&
        mv {params.trim_dir}/{wildcards.sample}_trimmed.fq.gz {output.trimmed}
        """


ruleorder: trimgalore_paired > trimgalore_single
