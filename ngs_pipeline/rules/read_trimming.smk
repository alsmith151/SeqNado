
def get_paired_fastq_inputs(wc):
    import glob
    fq_files = sorted(glob.glob(f"fastq/{wc.sample}_*.fastq.gz"))
    return {"fq1": fq_files[0], "fq2": fq_files[1]}

rule trimgalore_paired:
    # Trim reads using trimgalore
    input:
        unpack(get_paired_fastq_inputs),
    output:
        trimmed1=temp("trimmed/{sample}_1.fastq.gz"),
        trimmed2=temp("trimmed/{sample}_2.fastq.gz"),
    threads:
        4
    log:
        "logs/trimming/{sample}.log"
    shell: 
        """
           echo "trim_galore --cores {threads} --trim-n --paired --output_dir trimmed {input.fq1} {input.fq2}" > {log} &&
           trim_galore --cores {threads} --trim-n --paired --output_dir trimmed {input.fq1} {input.fq2} >> {log} 2>&1 &&
           mv trimmed/{wildcards.sample}_1_val_1.fq.gz {output.trimmed1} &&
           mv trimmed/{wildcards.sample}_2_val_2.fq.gz {output.trimmed2}
        """

rule trimgalore_single:
    # Trim reads using trimgalore
    input:
        fq1="fastq/{sample}.fastq.gz",
    output:
        trimmed1=temp("trimmed/{sample}.fastq.gz"),
    threads:
        4
    log:
        "logs/trimming/{sample}.log"
    shell: 
        """trim_galore --cores {threads} --trim-n --output_dir trimmed {input.fq1} > {log} 2>&1 &&
           mv trimmed/{wildcards.sample}_trimmed.fq.gz {output.trimmed1}
        """


ruleorder: trimgalore_paired > trimgalore_single