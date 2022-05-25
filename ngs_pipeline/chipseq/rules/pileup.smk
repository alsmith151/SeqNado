rule homer_make_tag_directory:
    input:
        bam = "aligned_and_filtered/{sample}.bam",
    output:
        homer_tag_directory = directory("tag_dirs/{sample}"),
    params:
        options = config["homer"]["maketagdirectory_options"] if config["homer"]["maketagdirectory_options"] else "",
    shell:
        """makeTagDirectory {output.homer_tag_directory} {input.bam} {params.options}"""


rule homer_make_bigwigs:
    input:
        homer_tag_directory = "tag_dirs/{sample}",
    output:
        homer_bigwig = "bigwigs/homer/{sample}.bigWig",
    log:
        log = "logs/homermakebigwigs_{sample}.log",
    params:
        genome_name = config["genome"]["name"],
        genome_chrom_sizes = config["genome"]["chrom_sizes"],
        options = config["homer"]["makebigwig_options"] if config["homer"]["makebigwig_options"] else "",
    run:
        
        cmd = f"""makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir bigwigs/homer/ {params.options} 2>&1 {log.log} &&
                 mv {output.homer_bigwig.replace(".bigWig", ".ucsc.bigWig")} {output.homer_bigwig}"""
        
        shell(cmd)

rule deeptools_make_bigwigs:
    input:
        bam = "aligned_and_filtered/{sample}.bam",
        bai = "aligned_and_filtered/{sample}.bam.bai",
    output:
        bigwig = "bigwigs/deeptools/{sample}.bigWig",
    params:
        threads = config["pipeline"]["n_cores"], 
        options = config["deeptools"]["bamcoverage_options"] if config["deeptools"]["bamcoverage_options"] else "",
    shell:
        """
        bamCoverage {params.options} -b {input.bam} -o {output.bigwig} -p {params.threads}
        """

