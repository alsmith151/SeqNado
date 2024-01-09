from seqnado import utils

rule fastq_screen:
    input:
        unpack(lambda wc: utils.translate_fq_files(wc, samples=FASTQ_SAMPLES, paired=False)),
    output:
        fq_screen=temp("seqnado_output/qc/fastq_screen/{sample}_{read}_screen.html"),
        fq_screen_png=temp("seqnado_output/qc/fastq_screen/{sample}_{read}_screen.png"),
        fq_screen_txt=temp("seqnado_output/qc/fastq_screen/{sample}_{read}_screen.txt"),
    params:
        outdir="seqnado_output/qc/fastq_screen",
        conf=config["spikein"]["fastq_screen_config"],
        tmpdir="seqnado_output/qc/fastqc_raw/{sample}_{read}",
        basename=lambda wc, output: utils.get_fq_filestem(wc, samples=FASTQ_SAMPLES),
    threads: config["bowtie2"]["threads"]
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * 2**attempt,
    log:
        "seqnado_output/logs/fastq_screen/{sample}_{read}.log",
    shell:"""
    fastq_screen --conf {params.conf} --threads {threads} --subset 10000 --aligner bowtie2 --threads 8 {input.fq}  --outdir {params.tmpdir} &&
        mv {params.tmpdir}/{params.basename}_screen.html {output.fq_screen} &&
        mv {params.tmpdir}/{params.basename}_screen.png {output.fq_screen_png} &&
        mv {params.tmpdir}/{params.basename}_screen.txt {output.fq_screen_txt} &&
        rm -r {params.tmpdir}
    """


rule multiqc_fastqscreen:
    input:
        expand(
            "seqnado_output/qc/fastq_screen/{sample}_{read}_screen.{ext}",
            sample=SAMPLE_NAMES,
            read=[1, 2],
            ext=["txt", "png", "html"]
        ),
    output:
        "seqnado_output/qc/full_fastqscreen_report.html",
    log:
        "seqnado_output/logs/multiqc_fastqscreen.log",
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * 2**attempt,
    shell:
        "multiqc -o seqnado_output/qc -n full_fastqscreen_report.html --force seqnado_output/qc/fastq_screen > {log} 2>&1"



rule align_paired_spikein:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
    params:
        index=config["genome"]["indicies"],
        options="--no-mixed --no-discordant",
    output:
        bam=temp("seqnado_output/aligned/spikein/raw/{sample}.bam"),
    threads: config["bowtie2"]["threads"]
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * 2**attempt,
    log:
        "seqnado_output/logs/align/{sample}.log",
    shell:"""
        bowtie2 -p {threads} {params.options} -x {params.index} -1 {input.fq1} -2 {input.fq2} 2> {log} |
        samtools view -bS - > {output.bam} &&
        samtools sort -@ {threads} -o {output.bam}_sorted {output.bam} >> {log} 2>&1 &&
        mv {output.bam}_sorted {output.bam}
        """


use rule sort_bam as sort_bam_spikein with:
    input:
        bam=rules.align_paired_spikein.output.bam,
    output:
        bam=temp("seqnado_output/aligned/spikein/sorted/{sample}.bam"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}_sort.log"
    
use rule index_bam as index_bam_spikein with:
    input:
        bam=rules.sort_bam_spikein.output.bam,
    output:
        bai=temp("seqnado_output/aligned/spikein/sorted/{sample}.bam.bai"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}_index.log"

rule filter_bam_spikein:
    input:
        bam=rules.sort_bam_spikein.output.bam,
    output:
        bam=temp("seqnado_output/aligned/spikein/filtered/{sample}.bam"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}_filter.log"
    shell:"""
    samtools view -b -F 3332 -q 30 -@ 8 {input.bam} > {output.bam} &&
    echo 'Filtered bam number of mapped reads:' > {log} 2>&1 &&
    samtools view -c {output.bam} >> {log} 2>&1
    """

use rule index_bam as index_bam_spikein_filtered with:
    input:
        bam=rules.filter_bam_spikein.output.bam,
    output:
        bai=temp("seqnado_output/aligned/spikein/filtered/{sample}.bam.bai"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}_filter_index.log"


rule split_bam:
    input:
        bam=rules.filter_bam_spikein.output.bam,
        bai=rules.index_bam_spikein_filtered.output.bai,
    output:
        ref_bam=temp("seqnado_output/aligned/spikein/{sample}_ref.bam"),
        exo_bam=temp("seqnado_output/aligned/spikein/{sample}_exo.bam"),
        stats="seqnado_output/aligned/spikein/{sample}_stats.tsv",
    params:
        genome_prefix=config["spikein"]["reference_genome"],
        exo_prefix=config["spikein"]["spikein_genome"],
        prefix="seqnado_output/aligned/spikein/{sample}",
        map_qual=30,
    log:"seqnado_output/logs/split_bam/{sample}.log"
    shell:"""
        samtools view -h {input.bam} | awk '{{if($0 ~ /^@/ || $3 ~ /^chr/) print}}' | samtools view -b -o {output.ref_bam} &&
        samtools view -h {input.bam} | awk '{{if($0 ~ /^@/ || $3 ~ /^{params.exo_prefix}/) print}}' | samtools view -b -o {output.exo_bam} &&
        echo -e "sample\treference_reads\tspikein_reads" > {output.stats} &&
        echo -e "{wildcards.sample}\t$(samtools view -c {output.ref_bam})\t$(samtools view -c {output.exo_bam})" >> {output.stats}
        """

rule move_ref_bam:
    input:
        bam=rules.split_bam.output.ref_bam,
    output:
        bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
    shell:"""
    mv {input.bam} {output.bam}
    """

rule calculate_normalisation_factors:
    input:
        expand(rules.split_bam.output.stats, sample=SAMPLE_NAMES),
    output:
        normalisation_table = "seqnado_output/normalisation_factors.tsv",
        normalisation_factors = "seqnado_output/normalisation_factors.json",
        scale_factors="seqnado_output/scaling_factors.json",
    log:
        "seqnado_output/logs/normalisation_factors.log"
    script:
        "../scripts/calculate_spikein_norm_orlando.py"