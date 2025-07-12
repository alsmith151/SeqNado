from seqnado.helpers import check_options, define_time_requested, define_memory_requested

if ASSAY == "ChIP":
    prefix = SAMPLE_NAMES_IP
else:
    prefix = SAMPLE_NAMES


rule make_dataset_regions:
    """Create a dataset from bigWig files using either a BED file."""
    input:
        bigwigs=expand(
            "seqnado_output/bigwigs/deeptools/{method}/{sample}.bigWig",
            method=get_scale_method(config) or "unscaled",
            sample=prefix,
        ),
    output:
        dataset="seqnado_output/dataset/dataset_regions.h5ad",
    params:
        bigwig_dir="seqnado_output/bigwigs/deeptools/unscaled/",
        chromosome_sizes=config['genome']['chromosome_sizes'],
        blacklist=check_options(config["blacklist"]),
        regions=check_options(config["dataset"]["regions_bed"]),
    threads: 1
    resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest",
    log: "seqnado_output/logs/make_dataset_regions.log",
    shell:"""
    quantnado-make-dataset \
    --bigwig-dir {params.bigwig_dir} \
    --output-file {output.dataset} \
    --chromsizes {params.chromosome_sizes} \
    --regions {params.regions} \
    --blacklist {params.blacklist} \
    --log-file {log}
    """

rule make_dataset_binsize:
    """Create a dataset from bigWig files using bin size."""
    input:
        bigwigs=expand(
            "seqnado_output/bigwigs/deeptools/{method}/{sample}.bigWig",
            method=get_scale_method(config) or "unscaled",
            sample=prefix,
        ),
    output:
        dataset="seqnado_output/dataset/dataset_bins.h5ad",
    params:
        bigwig_dir="seqnado_output/bigwigs/deeptools/unscaled/",
        chromosome_sizes=config['genome']['chromosome_sizes'],
        blacklist=check_options(config["blacklist"]),
        binsize=check_options(config["dataset"]["binsize"]),
    threads: 1
    resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest",
    log: "seqnado_output/logs/make_dataset_binsize.log",
    shell:"""
    quantnado-make-dataset \
    --bigwig-dir {params.bigwig_dir} \
    --output-file {output.dataset} \
    --chromsizes {params.chromosome_sizes} \
    --binsize {params.binsize} \
    --blacklist {params.blacklist} \
    --log-file {log}
    """
