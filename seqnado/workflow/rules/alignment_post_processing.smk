from seqnado.helpers import check_options, define_time_requested, define_memory_requested


rule sort_bam:
    input:
        bam="seqnado_output/aligned/raw/{sample}.bam",
    output:
        bam=temp("seqnado_output/aligned/sorted/{sample}.bam"),
        read_log=temp("seqnado_output/qc/alignment_post_process/{sample}_sort.tsv"),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=12, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    threads: config["samtools"]["threads"]
    log: "seqnado_output/logs/alignment_post_process/{sample}_sort.log",
    shell:"""
    samtools sort {input.bam} -@ {threads} -o {output.bam} -m 900M &&
    echo 'Step\tRead Count' > {output.read_log} &&
    echo -e "Raw counts\t$(samtools view -c {input.bam})" >> {output.read_log} &&
    echo -e "sort\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
    """


rule index_bam:
    input:
        bam="seqnado_output/aligned/sorted/{sample}.bam",
    output:
        bai=temp("seqnado_output/aligned/sorted/{sample}.bam.bai"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    shell:"samtools index -@ {threads} -b {input.bam}"


if config["remove_blacklist"] and os.path.exists(config.get("blacklist", "")):

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
            blacklist=check_options(config["blacklist"]),
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


if config["remove_pcr_duplicates_method"] == "picard":
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
        picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} --REMOVE_DUPLICATES true --CREATE_INDEX true {params.options} &&
        mv seqnado_output/aligned/duplicates_removed/{wildcards.sample}.bai {output.bai} &&
        echo -e "duplicate removal\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
        """

elif config["remove_pcr_duplicates_method"] == "samtools":
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


if config.get("shift_atac_reads"):
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
        rsbamtk shift -b {input.bam} -o {output.tmp} &&
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
    threads: config["samtools"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log: "seqnado_output/logs/alignment_post_process/{sample}_filter.log",
    params:
        options=check_options(config["samtools"]["filter_options"]),
    shell:"""
    samtools view -@ {threads} -h -b {input.bam} {params.options} > {output.bam} &&
    samtools index {output.bam} &&
    echo -e "filtering\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
    """

rule move_bam_to_final_location:
    input:
        bam="seqnado_output/aligned/filtered/{sample}.bam",
        bai="seqnado_output/aligned/filtered/{sample}.bam.bai",
    output:
        bam="seqnado_output/aligned/{sample,[A-Za-z\\d\\-_]+}.bam",
        bai="seqnado_output/aligned/{sample,[A-Za-z\\d\\-_]+}.bam.bai",
        read_log=temp("seqnado_output/qc/alignment_post_process/{sample}_final.tsv"),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log: "seqnado_output/logs/alignment_post_process/{sample}_final.log",
    shell:"""
    mv {input.bam} {output.bam} &&
    mv {input.bai} {output.bai} &&
    echo -e "Final reads\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
    """


rule bam_stats:
    input: 
        sort="seqnado_output/qc/alignment_post_process/{sample}_sort.tsv",
        blacklist="seqnado_output/qc/alignment_post_process/{sample}_blacklist.tsv",
        remove_duplicates="seqnado_output/qc/alignment_post_process/{sample}_remove_duplicates.tsv",
        atac_shift="seqnado_output/qc/alignment_post_process/{sample}_atac_shift.tsv",
        filtered="seqnado_output/qc/alignment_post_process/{sample}_filter.tsv",
        final="seqnado_output/qc/alignment_post_process/{sample}_final.tsv",
    output: temp("seqnado_output/qc/alignment_post_process/{sample}_alignment_stats.tsv")
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    shell: """
        cat {input.sort} {input.blacklist} {input.remove_duplicates} {input.atac_shift} {input.filtered} {input.final} > {output}
    """

rule prepare_stats_report:
    input:
        expand(
            "seqnado_output/qc/alignment_post_process/{sample}_alignment_stats.tsv",
            sample=SAMPLE_NAMES,
        ),
    output:
        "seqnado_output/qc/alignment_stats.tsv",
    log:
        "seqnado_output/logs/alignment_stats.log",
    script:
        "../scripts/alignment_stats.py"

def get_bam_files_for_merge(wildcards):
    from seqnado.design import NormGroups
    norm_groups = NormGroups.from_design(DESIGN, subset_column="merge")

    sample_names = norm_groups.get_grouped_samples(wildcards.group)
    bam_files = [
        f"seqnado_output/aligned/{sample}.bam" for sample in sample_names
    ]
    return bam_files


rule merge_bams:
    input:
        bams=get_bam_files_for_merge,
    output:
        temp("seqnado_output/aligned/merged/{group}.bam"),
    threads: config["samtools"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log: "seqnado_output/logs/merge_bam/{group}.log",
    shell:"""
    samtools merge {output} {input} -@ {threads}
    """


use rule index_bam as index_consensus_bam with:
    input:
        bam="seqnado_output/aligned/merged/{group}.bam",
    output:
        bai="seqnado_output/aligned/merged/{group}.bam.bai",
    threads: 8


localrules:
    move_bam_to_final_location,
