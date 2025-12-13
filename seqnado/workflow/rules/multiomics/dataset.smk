from seqnado.helpers import define_time_requested, define_memory_requested
from seqnado.config.multiomics import get_assay_bigwigs


bigwigs = get_assay_bigwigs(
    wildcards=None,
    ASSAYS=ASSAYS,
    rules=rules
)

SCALE_RESOURCES = 1

rule make_dataset_regions:
    """Create a dataset from bigWig files using either a BED file."""
    input:
        bigwigs=bigwigs,
    output:
        dataset=OUTPUT_DIR + "multiomics/dataset/dataset_regions.h5ad",
    params:
        bigwig_dir=OUTPUT_DIR + "multiomics/bigwigs/",
        chromosome_sizes=lambda wildcards: LOADED_CONFIGS[ASSAYS[0]]["genome"]["chromosome_sizes"],
        blacklist=lambda wildcards: LOADED_CONFIGS[ASSAYS[0]]["genome"]["blacklist"],
        regions=lambda wildcards: LOADED_CONFIGS[ASSAYS[0]]["assay_config"]["dataset_for_ml"]["regions_bed"],
    threads: 1
    resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest",
    log: OUTPUT_DIR + "multiomics/logs/make_dataset_regions.log",
    benchmark: OUTPUT_DIR + "multiomics/.benchmark/make_dataset_regions.tsv",
    message: "Making dataset from regions for machine learning"
    shell: """
    # symlink bigwig files to expected directory
    mkdir -p {params.bigwig_dir}
    for bw in {input.bigwigs}; do
        ln -s $(realpath $bw) {params.bigwig_dir}/$(basename $bw)
    done

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
        bigwigs=bigwigs,
    output:
        dataset=OUTPUT_DIR + "multiomics/dataset/dataset_bins.h5ad",
    params:
        bigwig_dir=OUTPUT_DIR + "multiomics/bigwigs/",
        chromosome_sizes=lambda wildcards: LOADED_CONFIGS[ASSAYS[0]]["genome"]["chromosome_sizes"],
        blacklist=lambda wildcards: LOADED_CONFIGS[ASSAYS[0]]["genome"]["blacklist"],
        binsize=lambda wildcards: LOADED_CONFIGS[ASSAYS[0]]["assay_config"]["dataset_for_ml"]["binsize"],
    threads: 1
    resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest",
    log: OUTPUT_DIR + "multiomics/logs/make_dataset_binsize.log",
    benchmark: OUTPUT_DIR + "multiomics/.benchmark/make_dataset_binsize.tsv",
    message: "Making dataset from binsize for machine learning"
    shell: """
    # symlink bigwig files to expected directory
    mkdir -p {params.bigwig_dir}
    for bw in {input.bigwigs}; do
        ln -s $(realpath $bw) {params.bigwig_dir}/$(basename $bw)
    done
    
    quantnado-make-dataset \
    --bigwig-dir {params.bigwig_dir} \
    --output-file {output.dataset} \
    --chromsizes {params.chromosome_sizes} \
    --binsize {params.binsize} \
    --blacklist {params.blacklist} \
    --log-file {log}
    """
