import seqnado.utils

def get_fq_files(wc):
    return {"fq1": FASTQ_SAMPLES.translation[f"{wc.sample}_1.fastq.gz"],
            "fq2": FASTQ_SAMPLES.translation[f"{wc.sample}_2.fastq.gz"]}

rule trimgalore_paired:
    # Trim reads using trimgalore
    input:
        unpack(get_fq_files)
    output:
        trimmed1=temp("seqnado_output/trimmed/{sample}_1.fastq.gz"),
        trimmed2=temp("seqnado_output/trimmed/{sample}_2.fastq.gz"),
    threads: 4
    params:
        options=seqnado.utils.check_options(config['trim_galore']['options'])
    log:
        "seqnado_output/logs/trimming/{sample}.log",
    shell:
        """
        trim_galore --cores {threads} {params.options} --paired --output_dir trimmed {input.fq1} {input.fq2} >> {log} 2>&1 &&
        mv trimmed/{wildcards.sample}_1_val_1.fq.gz {output.trimmed1} &&
        mv trimmed/{wildcards.sample}_2_val_2.fq.gz {output.trimmed2}
        """
