from seqnado.helpers import check_options, define_time_requested, define_memory_requested



rule fastq_screen_paired:
    input:
        fq="seqnado_output/fastqs/{sample}_{read}.fastq.gz",
    output:
        fq_screen=temp("seqnado_output/qc/fastq_screen/{sample}_{read}_screen.html"),
        fq_screen_png=temp("seqnado_output/qc/fastq_screen/{sample}_{read}_screen.png"),
        fq_screen_txt=temp("seqnado_output/qc/fastq_screen/{sample}_{read}_screen.txt"),
    params:
        outdir=temp("seqnado_output/qc/fastq_screen"),
        conf=config["fastq_screen_config"],
    threads: 4,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/fastq_screen/{sample}_{read}.log",
    shell:
        """ fastq_screen --conf {params.conf} --threads {threads} --subset 10000 --aligner bowtie2 --threads {threads} {input.fq}  --outdir {params.outdir} > {log} 2>&1 """


use rule fastq_screen_paired as fastq_screen_single with:
    input:
        fq="seqnado_output/fastqs/{sample}.fastq.gz",
    output:
        fq_screen=temp("seqnado_output/qc/fastq_screen/{sample}_screen.html"),
        fq_screen_png=temp("seqnado_output/qc/fastq_screen/{sample}_screen.png"),
        fq_screen_txt=temp("seqnado_output/qc/fastq_screen/{sample}_screen.txt"),
    log:
        "seqnado_output/logs/fastq_screen/{sample}.log",

ruleorder: fastq_screen_paired > fastq_screen_single 
