import seqnado.utils as utils


rule sort_bam:
    input:
        bam="aligned/raw/{sample}.bam",
    output:
        bam=temp("aligned/sorted/{sample}.bam"),
    threads: 8
    shell:
        """samtools sort {input.bam} -@ {threads} -o {output.bam} &&
           mv {input.bam}_sorted {input.bam}
        """


rule index_bam:
    input:
        bam="aligned/sorted/{sample}.bam",
    output:
        index=temp("aligned/sorted/{sample}.bam.bai"),
    threads: 1
    resources:
        mem_mb=1000
    shell:
        "samtools index -@ {threads} -b {input.bam}"


rule remove_blacklisted_regions:
    input:
        bam="aligned/sorted/{sample}.bam",
    output:
        bam=temp("aligned/blacklist_regions_removed/{sample}.bam"),
        bai=temp("aligned/blacklist_regions_removed/{sample}.bam.bai"),
    threads: 1
    params:
        blacklist=config["blacklist"],
    log:
        "logs/blacklist/{sample}.log",
    script:
        "../scripts/remove_blacklist.py"


rule remove_duplicates:
    input:
        bam="aligned/blacklist_regions_removed/{sample}.bam",
        index="aligned/blacklist_regions_removed/{sample}.bam.bai"
    output:
        bam=temp("aligned/duplicates_removed/{sample}.bam"),
        bai=temp("aligned/duplicates_removed/{sample}.bam.bai"),
    threads: 1
    log:
        "logs/duplicates/{sample}.log",
    script:
        "../scripts/remove_duplicates.py"

rule shift_atac_alignments:
    input:
        bam="aligned/blacklist_regions_removed/{sample}.bam",
    output:
        bam=temp("aligned/shifted_for_tn5_insertion/{sample}.bam"),
    params:
        options=None,
    threads: 1
    log:
        "logs/atac_shift/{sample}.log",
    script:
        "../scripts/shift_alignments.py"

rule move_bam_to_final_location:
    input:
        bam="aligned/shifted_for_tn5_insertion/{sample}.bam",
        bai="aligned/shifted_for_tn5_insertion/{sample}.bam.bai"
    output:
        bam="aligned/{sample}.bam",
        bai="aligned/{sample}.bam.bai",
    log:
        "logs/move_bam/{sample}.log",
    shell:
        """mv {input.bam} {output.bam} &&
           mv {input.bai} {output.bai} &&    
           echo "BAM moved to final location" > {log}
        """



# rule mark_filtering_complete:
#     input:
#         sentinel="flags/{sample}.blacklist.sentinel",
#         sentinel2="flags/{sample}.shifted.sentinel",
#     output:
#         sentinel=touch("flags/{sample}.filtering.complete.sentinel"),
#     log:
#         "logs/filtering/{sample}.log",
#     shell:
#         """echo "Filtering complete" > {log}"""

# localrules: mark_filtering_complete


# use rule index_bam as index_bam_filtered with:
#     input:
#         bam="aligned_and_filtered/{sample}.bam",
#         filtering_performed=rules.mark_filtering_complete.output.sentinel
#     output:
#         index="aligned_and_filtered/{sample}.bam.bai"
