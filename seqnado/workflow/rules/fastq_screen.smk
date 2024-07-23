
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
    threads: config.get("bowtie2")["threads"] if config.get("bowtie2") else 4,
    resources:
        mem=lambda wildcards, attempt: f"{16 * 2**attempt}GB",
    log:
        "seqnado_output/logs/fastq_screen/{sample}_{read}.log",
    shell:
        """ fastq_screen --conf {params.conf} --threads {threads} --subset 10000 --aligner bowtie2 --threads 8 {input.fq}  --outdir {params.outdir} > {log} 2>&1 """


use rule fastq_screen_paired as fastq_screen_single with:
    input:
        fq="seqnado_output/fastqs/{sample}.fastq.gz",
    output:
        fq_screen=temp("seqnado_output/qc/fastq_screen/{sample}_screen.html"),
        fq_screen_png=temp("seqnado_output/qc/fastq_screen/{sample}_screen.png"),
        fq_screen_txt=temp("seqnado_output/qc/fastq_screen/{sample}_screen.txt"),
    log:
        "seqnado_output/logs/fastq_screen/{sample}.log",


def get_fastqscreen_files(*args, **kwargs):
    """Return a list of fastq_screen files for a given sample name."""
    from pathlib import Path

    fastq_dir = Path("seqnado_output/fastqs/")
    fastqscreen_dir = Path("seqnado_output/qc/fastq_screen")
    fastqscreen_files = []
    for fq in fastq_dir.glob("*.fastq.gz"):
        base_name = fq.stem.split(".fastq")[0]
        for ext in ["html", "png", "txt"]:
            fastqscreen_file = fastqscreen_dir / f"{base_name}_screen.{ext}".replace(
                " ", ""
            )
            fastqscreen_files.append(str(fastqscreen_file))
    return fastqscreen_files


rule multiqc_fastqscreen:
    input:
        get_fastqscreen_files,
    output:
        "seqnado_output/qc/full_fastqscreen_report.html",
    log:
        "seqnado_output/logs/multiqc_fastqscreen.log",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        "multiqc -o seqnado_output/qc -n full_fastqscreen_report.html --force seqnado_output/qc/fastq_screen > {log} 2>&1"


ruleorder: fastq_screen_paired > fastq_screen_single > multiqc_fastqscreen
