import re
import seqnado.utils as utils
import pysam

# def is_sample_paired(sample):
#     if DESIGN.query("paired == True")["basename"].str.contains(sample).any():
#         return True


def is_bam_paired_end(wc, bam):

    if os.path.exists(bam):
        bam_ps = pysam.AlignmentFile(bam)
        head = bam_ps.head(1000)
        n_paired_reads = sum([aln.is_paired for aln in head])

        if n_paired_reads > 0:
            return True

    else:
        # TODO fix this for the new format
        # return is_sample_paired(wc.sample)
        return True


def filter_deeptools_bamcoverage_options(wc):

    bam = f"aligned_and_filtered/{wc.sample}.bam"
    options = (
        config["deeptools"]["bamcoverage"] if config["deeptools"]["bamcoverage"] else ""
    )

    if "-e" in options or "--extendReads" in options:
        if not is_bam_paired_end(wc, bam) and not re.search(
            "(-e \d+)|(--extendReads \d+)", options
        ):
            options = options.replace("--extendReads ", "").replace("-e ", "")

    return options


rule homer_make_tag_directory:
    input:
        bam="aligned_and_filtered/{sample}.bam",
        filtering="flags/{sample}.filtering.complete.sentinel",
    output:
        homer_tag_directory=directory("tag_dirs/{sample}"),
    params:
        options=config["homer"]["maketagdirectory"],
    log:
        "logs/homer/maketagdirectory_{sample}.log",
    shell:
        """makeTagDirectory {output.homer_tag_directory} {input.bam} {params.options} > {log} 2>&1"""


rule homer_make_bigwigs:
    input:
        homer_tag_directory="tag_dirs/{sample}",
        filtering="flags/{sample}.filtering.complete.sentinel",
    output:
        homer_bigwig="bigwigs/homer/{sample}.bigWig",
    log:
        "logs/homer/makebigwigs_{sample}.log",
    params:
        genome_name=config["genome"]["name"],
        genome_chrom_sizes=config["genome"]["chromosome_sizes"],
        options=config["homer"]["makebigwig"],
    run:
        cmd = f"""makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir bigwigs/homer/ {params.options} > {log} 2>&1 &&
                                 mv {output.homer_bigwig.replace(".bigWig", ".ucsc.bigWig")} {output.homer_bigwig}"""

        if workflow.use_singularity:
            cmd = utils.get_singularity_command(
                command=cmd,
                workflow=workflow,
            )

        shell(cmd)


rule deeptools_make_bigwigs:
    input:
        bam="aligned_and_filtered/{sample}.bam",
        bai="aligned_and_filtered/{sample}.bam.bai",
        filtering="flags/{sample}.filtering.complete.sentinel",
    output:
        bigwig="bigwigs/deeptools/{sample}.bigWig",
    params:
        options=lambda wc: filter_deeptools_bamcoverage_options(wc),
    threads: config["deeptools"]["threads"]
    log:
        "logs/pileups/deeptools/{sample}.log",
    shell:
        """
        bamCoverage {params.options} -b {input.bam} -o {output.bigwig} -p {threads} > {log} 2>&1
        """
