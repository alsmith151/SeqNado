import seqnado.utils as utils


rule sort_bam:
    input:
        bam="aligned/{sample}.bam",
    output:
        sentinel=touch("flags/{sample}.sorted.sentinel"),
    threads: 8
    shell:
        """samtools sort {input.bam} -@ {threads} -o {input.bam}_sorted &&
           mv {input.bam}_sorted {input.bam}
        """


rule index_bam:
    input:
        bam="aligned/{sample}.bam",
    output:
        index="aligned/{sample}.bam.bai",
    threads: 1
    resources:
        mem_mb=1000
    shell:
        "samtools index -@ {threads} -b {input.bam}"


rule remove_blacklisted_regions:
    input:
        bam="aligned_and_filtered/{sample}.bam",
    output:
        sentinel=touch("flags/{sample}.blacklist.sentinel"),
    threads: 1
    params:
        blacklist=config["blacklist"],
    log:
        "logs/blacklist/{sample}.log",
    script:
        "../scripts/remove_blacklist.py"


rule shift_atac_alignments:
    input:
        bam="aligned_and_filtered/{sample}.bam",
    output:
        sentinel=touch("flags/{sample}.shifted.sentinel"),
    params:
        options=None,
    threads: 1
    log:
        "logs/atac_shift/{sample}.log",
    script:
        "../scripts/shift_alignments.py"



rule mark_filtering_complete:
    input:
        sentinel="flags/{sample}.blacklist.sentinel",
        sentinel2="flags/{sample}.shifted.sentinel",
    output:
        sentinel=touch("flags/{sample}.filtering.complete.sentinel"),
    log:
        "logs/filtering/{sample}.log",
    shell:
        """echo "Filtering complete" > {log}"""

localrules: mark_filtering_complete


use rule index_bam as index_bam_filtered with:
    input:
        bam="aligned_and_filtered/{sample}.bam",
        filtering_performed=rules.mark_filtering_complete.output.sentinel
    output:
        index="aligned_and_filtered/{sample}.bam.bai"
