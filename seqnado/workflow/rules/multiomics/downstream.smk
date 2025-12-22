from seqnado.helpers import define_time_requested, define_memory_requested
from seqnado.outputs.multiomics import get_assay_bigwigs


# def get_bigwigs_for_dataset(wildcards):
#     """Get bigWig files at runtime after assay rules complete."""
#     return get_assay_bigwigs(wildcards, ASSAYS=ASSAYS, rules=rules)

rule gather_bigwigs:
    input:
        OUTPUT_DIR + "multiomics_summary.txt",
        bigwigs=MULTIOMICS_OUTPUT.bigwig_files,
    output:
        bw_dir = directory(OUTPUT_DIR + "multiomics/bigwigs/"),    
    threads: 1
    message: "Gathering bigWigs for multiomics downstream analyses"
    shell: """
    mkdir -p {output.bw_dir}
    for bw in {input.bigwigs}; do
        ln -sf $(realpath $bw) {output.bw_dir}/$(basename $bw)
    done
    """