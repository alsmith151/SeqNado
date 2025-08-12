
if CONFIG.shift_for_tn5_insertion:
    rule shift_atac_alignments:
        input:
            bam="seqnado_output/aligned/duplicates_removed/{sample}.bam",
            bai="seqnado_output/aligned/duplicates_removed/{sample}.bam.bai",
        output:
            bam=temp("seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam"),
            bai=temp(
                "seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam.bai"
            ),
            tmp=temp(
                "seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam.tmp"
            ),
            read_log=temp("seqnado_output/qc/alignment_post_process/{sample}_atac_shift.tsv"),
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        threads: 1
        log: "seqnado_output/logs/alignment_post_process/{sample}_atac_shift.log",
        shell:"""
        bamnado modify --input {input.bam} --output {output.tmp} --tn5-shift &&
        samtools sort {output.tmp} -@ {threads} -o {output.bam} &&
        samtools index {output.bam} &&
        echo -e "ATAC shift\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
        """

else:
    rule move_bam_to_temp_location:
        input:
            bam="seqnado_output/aligned/duplicates_removed/{sample}.bam",
            bai="seqnado_output/aligned/duplicates_removed/{sample}.bam.bai",
        output:
            bam=temp("seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam"),
            bai=temp(
                "seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam.bai"
            ),
            read_log=temp("seqnado_output/qc/alignment_post_process/{sample}_atac_shift.tsv"),
        threads: 1
        shell:"""
        mv {input.bam} {output.bam} &&
        mv {input.bam}.bai {output.bai} &&
        echo -e "ATAC shift\t$(samtools view -c {output.bam})" >> {output.read_log}
        """


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
    log: "seqnado_output/logs/alignment_post_process/{sample}_filter.log",
    params:
        options=str(CONFIG.third_party_tools.samtools.view.command_line_arguments),
    shell:"""
    samtools view -@ {threads} -h -b {input.bam} {params.options} > {output.bam} &&
    samtools index {output.bam} &&
    echo -e "filtering\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
    """