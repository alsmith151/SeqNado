if CONFIG.pcr_duplicates.tool == PCRDuplicateTool.PICARD:
    rule remove_duplicates_using_picard:
        input:
            bam="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam",
            bai="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam.bai",
        output:
            bam=temp("seqnado_output/aligned/duplicates_removed/{sample}.bam"),
            bai=temp("seqnado_output/aligned/duplicates_removed/{sample}.bam.bai"),
            metrics="seqnado_output/qc/library_complexity/{sample}.metrics",
            read_log=temp("seqnado_output/qc/alignment_post_process/{sample}_remove_duplicates.tsv"),
        threads: 8
        params:
            options=check_options(config["picard"]["options"]),
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=5, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        log: "seqnado_output/logs/alignment_post_process/{sample}_remove_duplicates.log",
        shell:"""
        picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} --REMOVE_DUPLICATES true --CREATE_INDEX true {params.options} 2> {log} &&
        mv seqnado_output/aligned/duplicates_removed/{wildcards.sample}.bai {output.bai} &&
        echo -e "duplicate removal\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
        """

elif CONFIG.pcr_duplicates.tool == PCRDuplicateTool.SAMTOOLS:
    rule remove_duplicates_using_samtools:
        input:
            bam="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam",
            bai="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam.bai",
        output:
            bam=temp("seqnado_output/aligned/duplicates_removed/{sample}.bam"),
            bai=temp("seqnado_output/aligned/duplicates_removed/{sample}.bam.bai"),
            read_log=temp("seqnado_output/qc/alignment_post_process/{sample}_remove_duplicates.tsv"),
        threads: config["samtools"]["threads"]
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=5, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        log: "seqnado_output/logs/alignment_post_process/{sample}_remove_duplicates.log",
        shell:"""
        samtools rmdup -@ {threads} {input.bam} {output.bam} &&
        samtools index {output.bam} &&
        echo -e "duplicate removal\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
        """
else:
    rule ignore_duplicates:
        input:
            bam="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam",
            bai="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam.bai",
        output:
            bam=temp("seqnado_output/aligned/duplicates_removed/{sample}.bam"),
            bai=temp("seqnado_output/aligned/duplicates_removed/{sample}.bam.bai"),
            read_log=temp("seqnado_output/qc/alignment_post_process/{sample}_remove_duplicates.tsv"),
        threads: 8
        resources:
            mem="500MB",
        log: "seqnado_output/logs/alignment_post_process/{sample}_remove_duplicates.log",
        shell: """
        mv {input.bam} {output.bam} &&
        mv {input.bai} {output.bai} &&
        echo -e "duplicate removal\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
        """

