import seqnado.utils as utils
PARTS=[str (x) for x in range(50)]
if config["split_fastq"] == "yes":
    if config["read_type"] == "paired":
        rule split_fq:
            input:
                unpack(lambda wc: seqnado.utils.translate_fq_files(wc, samples=FASTQ_SAMPLES, paired=True)),
            output:
                expand("seqnado_output/fastq_split/{{sample}}_{part}_{read}.fastq.gz", part=PARTS, read=["1", "2"]),
            params:
                split1=expand("-o seqnado_output/fastq_split/{{sample}}_{part}_1.fastq.gz", part=PARTS),
                split2=expand("-o seqnado_output/fastq_split/{{sample}}_{part}_2.fastq.gz", part=PARTS),
            resources:
                mem_mb=750,
            shell:"""
            fastqsplitter -i {input.fq1} {params.split1} &&
            fastqsplitter -i {input.fq2} {params.split2}
            """
        
        rule trimgalore_paired:
            input:
                split1="seqnado_output/fastq_split/{sample}_{part}_1.fastq.gz",
                split2="seqnado_output/fastq_split/{sample}_{part}_2.fastq.gz",
            output:
                trimmed1=temp("seqnado_output/trimmed/{sample}_{part}_1_trimmed.fq.gz"),
                trimmed2=temp("seqnado_output/trimmed/{sample}_{part}_2_trimmed.fq.gz"),
            threads: 4
            resources:
                mem_mb=750,
            params:
                options=utils.check_options(config['trim_galore']['options']),
                trim_dir="seqnado_output/trimmed"
            log:"seqnado_output/logs/trimming/{sample}_{part}.log",
            shell:"""
                trim_galore --cores {threads} {params.options} --basename {wildcards.sample}_{wildcards.part} --paired --output_dir {params.trim_dir} {input.split1} {input.split2} >> {log} 2>&1 &&
                mv {params.trim_dir}/{wildcards.sample}_{wildcards.part}_val_1.fq.gz {output.trimmed1} &&
                mv {params.trim_dir}/{wildcards.sample}_{wildcards.part}_val_2.fq.gz {output.trimmed2}
                """

        rule align_split:
            input:
                fq1="seqnado_output/trimmed/{sample}_{part}_1_trimmed.fq.gz",
                fq2="seqnado_output/trimmed/{sample}_{part}_2_trimmed.fq.gz",
            output:
                bam=temp("seqnado_output/aligned/split/{sample}_{part}.bam"),
            params:
                index=config["genome"]["indicies"],
                options=utils.check_options(config["bowtie2"]["options"]),
            threads: config["bowtie2"]["threads"]
            resources:
                mem_mb=4000 // int(config["bowtie2"]["threads"])
            log:"seqnado_output/logs/aligned/split/{sample}_part{part}.log",
            shell:"""
                bowtie2 -p {threads} -x {params.index} -1 {input.fq1} -2 {input.fq2} {params.options} 2> {log} |
                samtools view -bS - > {output.bam} &&
                samtools sort -@ {threads} -o {output.bam}_sorted {output.bam} >> {log} 2>&1 &&
                mv {output.bam}_sorted {output.bam}
                """

        rule merge_bams:
            input:
                expand("seqnado_output/aligned/split/{{sample}}_{part}.bam", part=PARTS),
            output:
                bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
            threads: 4
            log:"seqnado_output/logs/merge/{sample}.log",
            shell:"""
            samtools merge -o {output.bam} -@ {threads} -h {input} >> {log} 2>&1
            """
