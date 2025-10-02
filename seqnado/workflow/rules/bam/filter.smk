from seqnado.helpers import  define_time_requested, define_memory_requested

rule filter_bam:
    input:
        bam="seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam",
        bai="seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam.bai",
    output:
        bam="seqnado_output/aligned/filtered/{sample}.bam",
        bai="seqnado_output/aligned/filtered/{sample}.bam.bai",
        read_log=temp("seqnado_output/qc/alignment_post_process/{sample}_filter.tsv"),
    threads: CONFIG.third_party_tools.samtools.view.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: "seqnado_output/logs/alignment_post_process/{sample}_filter.log",
    params:
        options=str(CONFIG.third_party_tools.samtools.view.command_line_arguments),
    shell:"""
    samtools view -@ {threads} -h -b {input.bam} {params.options} > {output.bam} &&
    samtools index {output.bam} &&
    echo -e "filtering\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
    """
