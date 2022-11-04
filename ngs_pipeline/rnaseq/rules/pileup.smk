
# rule homer_make_bigwigs_rna:
#     input:
#         homer_tag_directory = "tag_dirs/{sample}",
#         filtering = "flags/{sample}.filtering.complete.sentinel",
#     output:
#         homer_bigwig = "bigwigs/homer/{sample}.bigWig",
#     log:
#         "logs/homer/makebigwigs_{sample}.log",
#     params:
#         genome_name = config["genome"]["name"],
#         genome_chrom_sizes = config["genome"]["chromosome_sizes"],
#         options = config["homer"]["makebigwig"],
#     run:
        
#         cmd = f"""makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir bigwigs/homer/ {params.options} > {log} 2>&1 &&
#                  mv {output.homer_bigwig.replace(".bigWig", ".ucsc.bigWig")} {output.homer_bigwig}"""
        
#         if workflow.use_singularity:
#             cmd = utils.get_singularity_command(command=cmd,
#                                                 workflow=workflow,)
        
#         shell(cmd)


# def filter_deeptools_bamcoverage_options(wc):
    
#     bam = f"aligned_and_filtered/{wc.sample}.bam"
#     options = config["deeptools"]["bamcoverage"] if config["deeptools"]["bamcoverage"] else ""

#     if "-e" in options or "--extendReads" in options:
#         if not is_bam_paired_end(wc, bam) and not re.search("(-e \d+)|(--extendReads \d+)", options):
#             options = options.replace("--extendReads ", "").replace("-e ", "")
    
#     return options

rule deeptools_make_bigwigs_rna_plus:
    input:
        bam = "aligned_and_filtered/{sample}.bam",
        bai = "aligned_and_filtered/{sample}.bam.bai",
        filtering = "flags/{sample}.filtering.complete.sentinel",
    output:
        bigwig = "bigwigs/deeptools/{sample}_plus.bigWig",
    params:
        options = config["deeptools"]["bamcoverage"],
    threads:
        config["deeptools"]["threads"],
    log:
        "logs/pileups/deeptools/{sample}.log"
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --filterRNAstrand forward -p {threads} > {log} 2>&1
        """

rule deeptools_make_bigwigs_rna_minus:
    input:
        bam = "aligned_and_filtered/{sample}.bam",
        bai = "aligned_and_filtered/{sample}.bam.bai",
        filtering = "flags/{sample}.filtering.complete.sentinel",
    output:
        bigwig = "bigwigs/deeptools/{sample}_minus.bigWig",
    params:
        options = config["deeptools"]["bamcoverage"],
    threads:
        config["deeptools"]["threads"],
    log:
        "logs/pileups/deeptools/{sample}.log"
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --filterRNAstrand reverse -p {threads} > {log} 2>&1
        """