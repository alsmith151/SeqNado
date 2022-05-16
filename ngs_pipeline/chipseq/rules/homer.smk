rule homer_make_tag_directory:
    input:
        bam = "aligned_and_postprocessed/{sample}.bam",
    output:
        homer_tag_directory = "tag_dirs/{sample}",
    params:
        options = config["homer"]["maketagdirectory_options"] if config["homer"]["maketagdirectory_options"] else "",
    shell:
        """makeTagDirectory {output.homer_tag_directory} {input.bam} {params.options}"""


rule homer_make_bigwigs:
    input:
        homer_tag_directory = "tag_dirs/{sample}",
    output:
        homer_bigwigs = "bigwigs/{sample}.bigwig" if not USE_DEEPTOOLS else "bigwigs/{sample}_homer.bigwig",
    params:
        genome_name = config["genome"]["name"],
        genome_chrom_sizes = config["genome"]["chrom_sizes"],
        options = config["homer"]["makebigwig_options"] if config["homer"]["makebigwig_options"] else "",
    shell:
        """makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir {output.homer_bigwigs} {params.options}"""
    