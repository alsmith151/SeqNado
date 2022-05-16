rule deeptools_filter_alignments:
    input:
        bam = "aligned/{sample}.bam",
    output:
        bam = "aligned_and_post-processed/{sample}.bam",
    params:
        threads = config["pipeline"]["n_cores"],
        deduplicate = "--ignoreDuplicates" if config["alignments"]["deduplicate"] else "",
        options = config["alignments"]["filter_options"] if config["alignments"]["filter_options"] else "",
    run:
        
        if params.deduplicate or params.options:

            shell("""
            alignmentSieve -b {input.bam} -o {output.bam} -p {params.threads} {params.deduplicate} {params.options} &&
            samtools sort -o {output.bam}.tmp {output.bam} -@ {params.threads} &&
            mv {output.bam}.tmp {output.bam} &&
            samtools index {output.bam} &&
            rm -f {output.bam}.tmp
            """)

        else:
            shell("""
                ln -s $(realpath {input.bam}) {output.bam} && 
                ln -s {input.bam}.bai {output.bam}.bai""")



rule deeptools_make_bigwigs:
    input:
        bam = "aligned_and_post-processed/{sample}.bam",
        bai = "aligned_and_post-processed/{sample}.bam.bai",
    output:
        bigwig = "bigwig/{sample}.bw" if not USE_HOMER else "bigwig/{sample}_deeptools.bigwig",
    params:
        threads = config["pipeline"]["n_cores"], 
        options = config["deeptools"]["bamcoverage_options"] if config["deeptools"]["bamcoverage_options"] else "",
    shell:
        """
        bamCoverage {params.options} -b {input.bam} -o {output.bigwig} -p {params.threads}
        """

