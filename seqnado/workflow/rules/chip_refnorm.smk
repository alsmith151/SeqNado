rule fastq_screen:
    input:
        unpack(lambda wc: seqnado.utils.translate_fq_files(wc, samples=FASTQ_SAMPLES, paired=False)),
    output:
        fq_screen="seqnado_output/qc/fastq_screen/{sample}_{read}_screen.html",
    params:
        outdir="seqnado_output/qc/fastq_screen",
        conf="/ceph/project/milne_group/shared/seqnado_reference/fastqscreen_reference/fastq_screen.conf",
    threads: config["bowtie2"]["threads"]
    resources:
        mem_mb=1500,
    log:
        "seqnado_output/logs/fastq_screen/{sample}_{read}.log",
    shell:"""
    fastq_screen --conf {params.conf} --threads {threads} --subset 10000 --aligner bowtie2 --threads 8 {input.fq}  --outdir {params.outdir}
    """


rule multiqc_fastqscreen:
    input:
        expand(
            "seqnado_output/qc/fastq_screen/{sample}_{read}_screen.html",
            sample=SAMPLE_NAMES,
            read=[1, 2],
        ),
    output:
        "seqnado_output/qc/full_fastqscreen_report.html",
    log:
        "seqnado_output/logs/multiqc_fastqscreen.log",
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * 2**attempt,
    shell:
        "multiqc -o seqnado_output/qc -n full_fastqscreen_report.html --force seqnado_output/qc/fastq_screen > {log} 2>&1"


use rule align_paired as align_paired_spikein with:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
    output:
        bam=temp("seqnado_output/aligned/spikein/raw/{sample}.bam"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}.log"


use rule sort_bam as sort_bam_spikein with:
    input:
        bam=rules.align_paired_spikein.output.bam,
    output:
        bam=temp("seqnado_output/aligned/spikein/sorted/{sample}.bam"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}.log"
    
use rule index_bam as index_bam_spikein with:
    input:
        bam=rules.sort_bam_spikein.output.bam,
    output:
        bai=temp("seqnado_output/aligned/spikein/sorted/{sample}.bam.bai"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}.index.log"

rule filter_bam_spikein:
    input:
        bam=rules.sort_bam_spikein.output.bam,
    output:
        bam=temp("seqnado_output/aligned/spikein/filtered/{sample}.bam"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}_filter.log"
    shell:"""
    samtools view -b -f 2 -F 260 -q 30  -@ 8 {input.bam} > {output.bam}
    """

use rule index_bam as index_bam_spikein_filtered with:
    input:
        bam=rules.filter_bam_spikein.output.bam,
    output:
        bai=temp("seqnado_output/aligned/spikein/filtered/{sample}.bam.bai"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}.index.log"


rule split_bam:
    input:
        bam=rules.filter_bam_spikein.output.bam,
        bai=rules.index_bam_spikein_filtered.output.bai,
    output:
        bam=temp("seqnado_output/aligned/spikein/{sample}_sample.bam"),
        clean_bam="seqnado_output/aligned/raw/{sample}.bam",
    params:
        genome_prefix=config["genome"]["sample_genome"],
        exo_prefix=config["genome"]["exo_genome"],
        prefix="seqnado_output/aligned/spikein/{sample}",
        map_qual=30,
    log:
        "seqnado_output/logs/split_bam/{sample}.log"
    script:
        """
        ../scripts/split_bam.py -i {input.bam} -o {params.prefix} -g {params.genome_prefix} -p {params.exo_prefix} -q {params.map_qual} > {log} 2>&1 &&
        mv {output.bam} > {output.clean_bam}
        """
