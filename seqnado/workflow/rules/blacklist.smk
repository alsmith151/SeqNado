
if CONFIG.remove_blacklist:

    rule remove_blacklisted_regions:
        input:
            bam="seqnado_output/aligned/sorted/{sample}.bam",
            bai=rules.index_bam.output.bai,
        output:
            bam=temp("seqnado_output/aligned/blacklist_regions_removed/{sample}.bam"),
            bai=temp(
                "seqnado_output/aligned/blacklist_regions_removed/{sample}.bam.bai"
            ),
            read_log=temp("seqnado_output/qc/alignment_post_process/{sample}_blacklist.tsv"),
        threads: 1
        params:
            blacklist=CONFIG.genome.blacklist,
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=5, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        log: "seqnado_output/logs/alignment_post_process/{sample}_blacklist.log",
        shell:"""
        bedtools intersect -v -b {params.blacklist} -a {input.bam} > {output.bam} &&
        samtools index -b {output.bam} -o {output.bai} &&
        echo -e "blacklisted regions removal\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
        """

else:

    rule ignore_blacklisted_regions:
        input:
            bam="seqnado_output/aligned/sorted/{sample}.bam",
            bai=rules.index_bam.output.bai,
        output:
            bam=temp("seqnado_output/aligned/blacklist_regions_removed/{sample}.bam"),
            bai=temp(
                "seqnado_output/aligned/blacklist_regions_removed/{sample}.bam.bai"
            ),
            read_log=temp("seqnado_output/qc/alignment_post_process/{sample}_blacklist.tsv"),
        threads: 1
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        log: "seqnado_output/logs/alignment_post_process/{sample}_blacklist.log",
        shell:"""
        mv {input.bam} {output.bam} &&
        mv {input.bai} {output.bai} &&
        echo -e "blacklisted regions removal\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
        """