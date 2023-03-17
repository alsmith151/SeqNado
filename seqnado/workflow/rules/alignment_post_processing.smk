import seqnado.utils as utils


rule sort_bam:
    input:
        bam="seqnado_output/aligned/raw/{sample}.bam",
    output:
        bam=temp("seqnado_output/aligned/sorted/{sample}.bam"),
    threads: 8
    shell:
        """samtools sort {input.bam} -@ {threads} -o {output.bam}
        """


rule index_bam:
    input:
        bam="seqnado_output/aligned/sorted/{sample}.bam",
    output:
        bai=temp("seqnado_output/aligned/sorted/{sample}.bam.bai"),
    threads: 1
    resources:
        mem_mb=1000
    shell:
        "samtools index -@ {threads} -b {input.bam}"


rule remove_blacklisted_regions:
    input:
        bam="seqnado_output/aligned/sorted/{sample}.bam",
        bai="seqnado_output/aligned/sorted/{sample}.bam.bai"
    output:
        bam=temp("seqnado_output/aligned/blacklist_regions_removed/{sample}.bam"),
        bai=temp("seqnado_output/aligned/blacklist_regions_removed/{sample}.bam.bai"),
    threads: 1
    params:
        blacklist=config["blacklist"],
    log:
        "seqnado_output/logs/blacklist/{sample}.log",
    script:
        "../scripts/remove_blacklist.py"


rule remove_duplicates:
    input:
        bam="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam",
        bai="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam.bai"
    output:
        bam=temp("seqnado_output/aligned/duplicates_removed/{sample}.bam"),
        bai=temp("seqnado_output/aligned/duplicates_removed/{sample}.bam.bai"),
    threads: 1
    log:
        "seqnado_output/logs/duplicates/{sample}.log",
    script:
        "../scripts/remove_duplicates.py"

rule shift_atac_alignments:
    input:
        bam="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam",
        bai="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam.bai"
    output:
        bam=temp("seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam"),
        bai=temp("seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam.bai"),
    params:
        options=None,
    threads: 1
    log:
        "seqnado_output/logs/atac_shift/{sample}.log",
    script:
        "../scripts/shift_alignments.py"

rule move_bam_to_final_location:
    input:
        bam="seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam",
        bai="seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam.bai"
    output:
        bam="seqnado_output/aligned/{sample,[A-Za-z0-9_\-]+}.bam",
        bai="seqnado_output/aligned/{sample,[A-Za-z0-9_\-]+}.bam.bai",
    log:
        "seqnado_output/logs/move_bam/{sample}.log",
    shell:
        """mv {input.bam} {output.bam} &&
           mv {input.bai} {output.bai} &&    
           echo "BAM moved to final location" > {log}
        """

localrules: move_bam_to_final_location
