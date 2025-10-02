rule save_design:
    output:
        "seqnado_output/design.csv",
    container:
        None
    run:
        DESIGN.to_dataframe().to_csv("seqnado_output/design.csv", index=False)


rule make_genomic_bins:
    input:
        chrom_sizes=CONFIG.genome.chromosome_sizes,
    params:
        bin_size=CONFIG.genome.bin_size,
    output:
        bed="seqnado_output/resources/genomic_bins.bed",
    log:
        "seqnado_output/logs/genomic_bins.log",
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        """
        cooler makebins {input.chrom_sizes} {params.bin_size} -o {output.bed} > {log} 2>&1
        """


rule bed_to_saf:
    input:
        bed="seqnado_output/resources/genomic_bins.bed",
    output:
        saf="seqnado_output/resources/genomic_bins.saf",
    log:
        "seqnado_output/logs/genomic_bins_saf.log",
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        """
        awk 'BEGIN{{OFS="\t"}} {{print $4, $1, $2, $3, "."}}' {input.bed} > {output.saf} 2> {log}
        """
