
def get_paired_fastq_inputs(wc):
    import glob
    fq_files = glob.glob(f"fastq/{wc.sample}_*.fastq.gz")
    return {"fq1": fq_files[0], "fq2": fq_files[1]}

rule trimgalore_paired:
    # Trim reads using trimgalore
    input:
        unpack(get_paired_fastq_inputs),
    output:
        trimmed1="trimmed/{sample}_1.fastq.gz",
        trimmed2="trimmed/{sample}_2.fastq.gz",
    threads:
        config["pipeline"]["n_cores"],
    log:
        "logs/trimgalore_paired_{sample}.log"
    shell: 
        """trim_galore --cores {threads} --trim-n --paired --output_dir trimmed {input.fq1} {input.fq2} 2> {log} &&
           mv trimmed/{wildcards.sample}_1_val_1.fq.gz {output.trimmed1} &&
           mv trimmed/{wildcards.sample}_2_val_2.fq.gz {output.trimmed2}
        """

rule trimgalore_single:
    # Trim reads using trimgalore
    input:
        fq1="fastq/{sample}.fastq.gz",
    output:
        trimmed1="trimmed/{sample}.fastq.gz",
    params:
        threads = config["pipeline"]["n_cores"],
    log:
        "logs/trimgalore_{sample}.log"
    shell: 
        """trim_galore --cores {params.threads} --trim-n --output_dir trimmed {input.fq1} 2> {log} &&
           mv trimmed/{wildcards.sample}_trimmed.fq.gz {output.trimmed1}
        """


ruleorder: trimgalore_paired > trimgalore_single