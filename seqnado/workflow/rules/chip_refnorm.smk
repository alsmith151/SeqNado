from seqnado import utils

rule fastq_screen:
    input:
        unpack(lambda wc: utils.translate_fq_files(wc, samples=FASTQ_SAMPLES, paired=False)),
    output:
        fq_screen="seqnado_output/qc/fastq_screen/{sample}_{read}_screen.html",
        fq_screen_png="seqnado_output/qc/fastq_screen/{sample}_{read}_screen.png",
        fq_screen_txt="seqnado_output/qc/fastq_screen/{sample}_{read}_screen.txt",
    params:
        outdir="seqnado_output/qc/fastq_screen",
        conf=config["genome"]["fastq_screen_config"],
        tmpdir="seqnado_output/qc/fastqc_raw/{sample}_{read}",
        basename=lambda wc, output: utils.get_fq_filestem(wc, samples=FASTQ_SAMPLES),
    threads: config["bowtie2"]["threads"]
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * 2**attempt,
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
            "seqnado_output/qc/fastq_screen/{sample}_{read}_screen.txt",
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
        ref_bam="seqnado_output/aligned/spikein/{sample}_ref.bam",
        exo_bam="seqnado_output/aligned/spikein/{sample}_exo.bam",
        bam="seqnado_output/aligned/raw/{sample}_exo.bam",
    params:
        genome_prefix=config["genome"]["sample_genome"],
        exo_prefix=config["genome"]["exo_genome"],
        prefix="seqnado_output/aligned/spikein/{sample}",
        map_qual=30,
    log:
        "seqnado_output/logs/split_bam/{sample}.log"
    shell:"""
        samtools view -h {input.bam} | awk '{{if($0 ~ /^@/ || $3 ~ /^chr/) print}}' | samtools view -b -o {output.ref_bam} &&
        samtools view -h {input.bam} | awk '{{if($0 ~ /^@/ || $3 ~ /^{params.exo_prefix}/) print}}' | samtools view -b -o {output.exo_bam} &&
        cp {output.ref_bam} {output.bam} &&
        echo 'reference_reads:' > {log} 2>&1 &&
        samtools view -F 0x04 -c {output.ref_bam} >> {log} 2>&1 &&
        echo 'spikein_reads:' >> {log} 2>&1 &&
        samtools view -F 0x04 -c {output.exo_bam} >> {log} 2>&1
        """

use rule index_bam as index_bam_split with:
    input:
        bam=rules.split_bam.output.ref_bam,
    output:
        bai=temp("seqnado_output/aligned/spikein/{sample}_ref.bam.bai"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}.index.log"

use rule samtools_stats as samtools_stats_reference with:
    input:
        bam=rules.split_bam.output.ref_bam,
    output:
        stats="seqnado_output/aligned/spikein/{sample}_ref.txt",

use rule samtools_stats as samtools_stats_exo with:
    input:
        bam=rules.split_bam.output.exo_bam,
    output:
        stats="seqnado_output/aligned/spikein/{sample}_exo.txt",

# rule calculate_normalisation_factors:
#     input:
#         bam_ref = expand(bam=rules.split_bam.output.ref_bam),
#         bam_spikein = expand(bam=rules.split_bam.output.exo_bam),
#         bam_spikein_index = expand(bam=rules.index_bam_split.output.bai),
#         design = rules.save_design.output[0],
#     output:
#         normalisation_table = "seqnado_output/normalisation_factors.tsv",
#         # normalisation_factors = "seqnado_output/normalisation_factors.json"
#     container:
#         None
#     log:
#         "seqnado_output/logs/normalisation_factors.log"
#     script:
#         "../scripts/calculate_spikein_norm_factors.py"
