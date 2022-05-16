rule samtools_stats_raw:
    input:
        bam = "aligned/{sample}.bam"
    output:
        stats = "statistics/alignment/{sample}.txt"
    shell:
       """samtools stats {input.bam} > {output.stats}"""

rule multiqc_raw:
    input:
        stats = expand("statistics/alignment/{sample}.txt", sample=SAMPLES)
    output:
        report = "statistics/alignmentqc_report.html"
    shell:
        "multiqc statistics/alignment/ -o statistics -n alignmentqc_report.html"
