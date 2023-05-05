import seqnado.utils

rule trimgalore_paired:
    # Trim reads using trimgalore
    input:
        unpack(lambda wc: seqnado.utils.translate_fq_files(wc, samples=FASTQ_SAMPLES, paired=True))
    output:
        trimmed1=temp("seqnado_output/trimmed/{sample}_1.fastq.gz"),
        trimmed2=temp("seqnado_output/trimmed/{sample}_2.fastq.gz"),
    threads: 4
    resources:
        mem_mb=750,
        time="02:00:00",
    params:
        options=seqnado.utils.check_options(config['trim_galore']['options']),
        trim_dir="seqnado_output/trimmed"
    log:
        "seqnado_output/logs/trimming/{sample}.log",
    shell:
        """
        trim_galore --cores {threads} {params.options} --basename {wildcards.sample} --paired --output_dir {params.trim_dir} {input.fq1} {input.fq2} >> {log} 2>&1 &&
        mv {params.trim_dir}/{wildcards.sample}_val_1.fq.gz {output.trimmed1} &&
        mv {params.trim_dir}/{wildcards.sample}_val_2.fq.gz {output.trimmed2}
        """
