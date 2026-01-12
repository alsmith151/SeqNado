reference_genome = config["spikein_options"]["reference_genome"]
spikein_genome = config["spikein_options"]["spikein_genome"]


rule rename_chr_cat_fasta:
    input:
        reference_fasta="/project/milne_group/cchahrou/reference/genome/{reference_genome}/{reference_genome}.fa.gz",
        spikein_fasta="/project/milne_group/cchahrou/reference/genome/{spikein_genome}/{spikein_genome}.fa.gz",
    output:
        bt2_index="data/index_bt2/{reference_genome}_{spikein_genome}.1.bt2",
    params:
        cat_genome="data/{reference_genome}_{spikein_genome}.fa",
        prefix="data/index_bt2/{reference_genome}_{spikein_genome}",
        spikein_genome=spikein_genome,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    threads: 4
    log: OUTPUT_DIR + "/logs/genomes/rename_chr_cat_fasta.log",
    benchmark: OUTPUT_DIR + "/.benchmark/genomes/rename_chr_cat_fasta.tsv",
    message: "Renaming chromosomes and concatenating FASTA files",
    shell: """
    gunzip {input.spikein_fasta} | sed -e 's/chr/{params.spikein_genome}/' | cat -e {input.reference_fasta} > {params.cat_genome} &&
    bowtie2-build {params.cat_genome} {params.prefix}
    """


rule STAR_index:
    input:
        reference_fasta="/project/milne_group/cchahrou/reference/genome/{reference_genome}/{reference_genome}.fa.gz",
        spikein_fasta="/project/milne_group/cchahrou/reference/genome/{spikein_genome}/{spikein_genome}.fa.gz",
        gtf_in=config["gtf"],
    output:
        star_index="data/index_STAR/{reference_genome}_{spikein_genome}",
    params:
        prefix="data/index_STAR/{reference_genome}_{spikein_genome}",
        spikein_genome=spikein_genome,
    threads: 6
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/genomes/STAR_index.log",
    benchmark: OUTPUT_DIR + "/.benchmark/genomes/STAR_index.tsv",
    message: "Generating STAR index with spike-in sequences",
    shell:
        """
    gunzip {input.spikein_fasta} | sed -e 's/chr/{params.spikein_genome}/' | cat -e {input.reference_fasta} > {params.cat_genome} &&
    STAR --version &&
    STAR --runThreadN 6 \
    --runMode genomeGenerate \
    --genomeDir {params.prefix} \
    --genomeFastaFiles {input.reference_in} \
    --sjdbGTFfile {input.gtf_in} \
    --sjdbOverhang 99
    """
