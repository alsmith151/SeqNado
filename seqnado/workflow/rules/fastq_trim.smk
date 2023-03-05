import seqnado.utils

rule trimgalore_paired:
    # Trim reads using trimgalore
    input:
        fq1="fastq/{sample}_1.fastq.gz",
        fq2="fastq/{sample}_2.fastq.gz",
    output:
        trimmed1=temp("trimmed/{sample}_1.fastq.gz"),
        trimmed2=temp("trimmed/{sample}_2.fastq.gz"),
    threads: 4
    params:
        options=seqnado.utils.check_options(config['trim_galore']['options'])
    log:
        "logs/trimming/{sample}.log",
    shell:
        """
        trim_galore --cores {threads} {params.options} --paired --output_dir trimmed {input.fq1} {input.fq2} >> {log} 2>&1 &&
        mv trimmed/{wildcards.sample}_1_val_1.fq.gz {output.trimmed1} &&
        mv trimmed/{wildcards.sample}_2_val_2.fq.gz {output.trimmed2}
        """


rule trimgalore_single:
    # Trim reads using trimgalore
    input:
        fq1="fastq/{sample}.fastq.gz",
    output:
        trimmed1=temp("trimmed/{sample}.fastq.gz"),
    threads: 4
    params:
        options=seqnado.utils.check_options(config['trim_galore']['options'])
    log:
        "logs/trimming/{sample}.log",
    shell:
        """trim_galore --cores {threads} {params.options} --output_dir trimmed {input.fq1} > {log} 2>&1 &&
           mv trimmed/{wildcards.sample}_trimmed.fq.gz {output.trimmed1}
        """


ruleorder: trimgalore_paired > trimgalore_single
