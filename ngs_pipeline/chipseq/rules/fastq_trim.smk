rule trimgalore_paired:
    # Trim reads using trimgalore
    input:
        fq1="fastq/{sample}_R1.fastq.gz",
        fq2="fastq/{sample}_R2.fastq.gz",
    output:
        trimmed1="trimmed/{sample}_R1.fastq.gz",
        trimmed2="trimmed/{sample}_R2.fastq.gz",
    params:
        threads = config["pipeline"]["n_cores"],
    shell: 
        """trim_galore --cores {params.threads} --trim-n --paired --output_dir trimmed {input.fq1} {input.fq2} &&
           mv trimmed/{wildcards.sample}_1_val_1.fq.gz {output.trimmed1} &&
           mv trimmed/{wildcards.sample}_2_val_2.fq.gz {output.trimmed2}
        """

rule trimgalore_single:
    # Trim reads using trimgalore
    input:
        fq1="fastq/{sample}_R0.fastq.gz",
    output:
        trimmed1="trimmed/{sample}_R0.fastq.gz",
    params:
        threads = config["pipeline"]["n_cores"],
    shell: 
        """trim_galore --cores {params.threads} --trim-n --output_dir trimmed {input.fq1} &&
           mv trimmed/{wildcards.sample}_R0_trimmed.fq.gz {output.trimmed1}
        """


ruleorder: trimgalore_paired > trimgalore_single