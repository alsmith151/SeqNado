import seqnado.utils as utils


rule sort_bam:
    input:
        bam="aligned/raw/{sample}.bam",
    output:
        bam=temp("aligned/sorted/{sample}.bam"),
    threads: 8
    shell:
        """samtools sort {input.bam} -@ {threads} -o {output.bam}
        """


rule index_bam:
    input:
        bam="aligned/sorted/{sample}.bam",
    output:
        bai=temp("aligned/sorted/{sample}.bam.bai"),
    threads: 1
    resources:
        mem_mb=1000
    shell:
        "samtools index -@ {threads} -b {input.bam}"


rule remove_blacklisted_regions:
    input:
        bam="aligned/sorted/{sample}.bam",
        bai="aligned/sorted/{sample}.bam.bai"
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
        bai="aligned/blacklist_regions_removed/{sample}.bam.bai"
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
        bai="aligned/blacklist_regions_removed/{sample}.bam.bai"
    output:
        bam=temp("aligned/shifted_for_tn5_insertion/{sample}.bam"),
        bai=temp("aligned/shifted_for_tn5_insertion/{sample}.bam.bai"),
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
        bam="aligned/{sample,[A-Za-z0-9_\-]+}.bam",
        bai="aligned/{sample,[A-Za-z0-9_\-]+}.bam.bai",
    log:
        "logs/move_bam/{sample}.log",
    shell:
        """mv {input.bam} {output.bam} &&
           mv {input.bai} {output.bai} &&    
           echo "BAM moved to final location" > {log}
        """

localrules: move_bam_to_final_location
