from seqnado.helpers import  define_time_requested, define_memory_requested

rule sort_bam:
    input:
        bam=OUTPUT_DIR + "/aligned/raw/{sample}.bam",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/sorted/{sample}.bam"),
        read_log=temp(OUTPUT_DIR + "/qc/alignment_post_process/{sample}_sort.tsv"),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    threads: CONFIG.third_party_tools.samtools.sort.threads
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_sort.log",
    shell: """
        samtools sort {input.bam} -@ {threads} -o {output.bam} -m 900M
        echo 'Step\tRead Count' > {output.read_log}
        echo -e "Raw counts\t$(samtools view -c {input.bam})" >> {output.read_log}
        echo -e "sort\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
    """


rule index_bam:
    input:
        bam=OUTPUT_DIR + "/aligned/sorted/{sample}.bam",
    output:
        bai=temp(OUTPUT_DIR + "/aligned/sorted/{sample}.bam.bai"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:"samtools index -@ {threads} -b {input.bam}"


rule move_bam_to_final_location:
    input:
        bam=OUTPUT_DIR + "/aligned/filtered/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/filtered/{sample}.bam.bai",
    output:
        bam=OUTPUT_DIR + "/aligned/{sample,[A-Za-z\\d\-_]+}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample,[A-Za-z\\d\-_]+}.bam.bai",
        read_log=temp(OUTPUT_DIR + "/qc/alignment_post_process/{sample}_final.tsv"),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_final.log",
    shell:"""
    mv {input.bam} {output.bam} &&
    mv {input.bai} {output.bai} &&
    echo -e "Final reads\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
    """


localrules:
    move_bam_to_final_location