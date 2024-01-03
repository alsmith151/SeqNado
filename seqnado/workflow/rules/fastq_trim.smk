import seqnado.utils


def get_linked_fastq_files(wc):
    outdir = checkpoints.link_fastqs.get(**wc).output[0]
    fq_files = list(pathlib.Path(outdir).glob(f"{wc['sample']}*.fastq.gz"))
    return fq_files


checkpoint link_fastqs:
    input:
        fastqs=DESIGN.fastq_paths,
    output:
        directory("seqnado_output/fastqs"),
    run:
        for assay_name, assay in DESIGN.assays.items():
            if assay.r2:
                r1_path_new = pathlib.Path(
                    f"seqnado_output/fastqs/{assay_name}_1.fastq.gz"
                )
                r2_path_new = pathlib.Path(
                    f"seqnado_output/fastqs/{assay_name}_2.fastq.gz"
                )
                r1_path_new.symlink_to(assay.r1)
                r2_path_new.symlink_to(assay.r2)
            else:
                r1_path_new = pathlib.Path(
                    f"seqnado_output/fastqs/{assay_name}.fastq.gz"
                )
                r1_path_new.symlink_to(assay.r1)


rule confirm_fastqs_renamed:
    input:
        fastqs=get_linked_fastq_files,
    output:
        touch("seqnado_output/fastqs/.sentinel"),


rule trimgalore_paired:
    # Trim reads using trimgalore
    input:
        fq1="seqnado_output/fastqs/{sample}_1.fastq.gz",
        fq2="seqnado_output/fastqs/{sample}_2.fastq.gz",
        renamed_sentinel="seqnado_output/fastqs/.sentinel",
    output:
        trimmed1=temp("seqnado_output/trimmed/{sample}_1.fastq.gz"),
        trimmed2=temp("seqnado_output/trimmed/{sample}_2.fastq.gz"),
    threads: 4
    resources:
        mem_mb=2000,
        time="02:00:00",
    params:
        options=seqnado.utils.check_options(config["trim_galore"]["options"]),
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
        fq="seqnado_output/fastq/{sample}.fastq.gz",
        renamed_sentinel="seqnado_output/fastqs/.sentinel",
    output:
        trimmed=temp("seqnado_output/trimmed/{sample}.fastq.gz"),
    threads: 4
    resources:
        mem_mb=2000,
        time="02:00:00",
    params:
        options=seqnado.utils.check_options(config["trim_galore"]["options"]),
        trim_dir="seqnado_output/trimmed",
    log:
        "seqnado_output/logs/trimming/{sample}.log",
    shell:
        """
        trim_galore --cores {threads} {params.options} --basename {wildcards.sample} --output_dir {params.trim_dir} {input.fq} >> {log} 2>&1 &&
        mv {params.trim_dir}/{wildcards.sample}_trimmed.fq.gz {output.trimmed}
        """


ruleorder: trimgalore_paired > trimgalore_single
