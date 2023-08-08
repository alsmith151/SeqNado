
rule rename_chr_cat_fasta: # renames the chromosome names in the custom fasta file to not match the reference genome
    input: 
        custom_fasta = "/project/milne_group/cchahrou/reference/genome/dm6/dm6.fa.gz",
        reference_fasta = "/project/milne_group/cchahrou/reference/genome/hg38/hg38.fa",
    output: 
        ref_custom_genome = "data/ref_custom_genome.fasta",
    params:
        chr_prefix = "custom_chr",
    shell: """
        gunzip {input.custom_fasta} | sed -e 's/chr/{params.chr_prefix}/' | cat -e {input.reference_fasta} > {output.ref_custom_genome}
        """

rule bowtie2_index:
# indexes the custom genome for bowtie2
    input: 
        reference_in = rules.rename_chr_cat_fasta.output.ref_custom_genome,
    output:
        bt2_index = "data/index_bt2/ref_custom_genome.1.bt2",
    params:
        prefix = "data/index_bt2/ref_custom_genome",
    threads: config["bowtie2"]["threads"],
    shell: "bowtie2-build {input.reference_in} {params.prefix}"

# rule STAR_index:
# # indexes the custom genome for STAR
#     input: 
#         reference_in = rules.rename_chr_cat_fasta.output.ref_custom_genome,
#         gtf_in = config["gtf"],
#     output:
#         star_index = "data/index_STAR/ref_custom_genome",
#     params:
#         prefix = "data/index_STAR/ref_custom_genome",
#     threads: 6,
#     shell: """STAR --version && 
#     STAR --runThreadN 6 \
#     --runMode genomeGenerate \
#     --genomeDir {params.prefix} \
#     --genomeFastaFiles {input.reference_in} \
#     --sjdbGTFfile {input.gtf_in} \
#     --sjdbOverhang 99"""