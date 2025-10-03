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

rule get_fasta:
    input:
        peaks="seqnado_output/peaks/{method}/{sample}.bed",
    output:
        fasta="seqnado_output/motifs/fasta/{sample}.fa",
        bed=temp("motifs/fasta/{sample}_clean.bed"),
    params:
        genome=CONFIG.genome.fasta,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/motifs/fasta/{sample}.log",
    shell:
        """
    cat {input.peaks} | cut -f 1-3 > {output.bed} &&
    bedtools getfasta -fullHeader -fi {params.genome} -bed {output.bed} -fo {output.fasta}
    """

rule validate_peaks:
    input:
        peaks=OUTPUT.peaks,
    output:
        sentinel="seqnado_output/peaks/.validated",
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/validate_peaks.log",
    run:
        from loguru import logger

        with logger.catch():
            for peak_file in input.peaks:
                with open(peak_file, "r+") as p:
                    peak_entries = p.readlines()
                    n_peak_lines = sum(1 for line in peak_entries if not line.startswith("#"))
                    if n_peak_lines == 0: # empty peak file, write a dummy peak
                        p.write("chr21\t1\t2\n")

        with open(output.sentinel, "w") as s:
            s.write("validated")


rule bed_to_bigbed:
    input:
        bed="seqnado_output/peaks/{directory}/{sample}.bed",
        sentinel="seqnado_output/peaks/.validated",
    output:
        bigbed="seqnado_output/peaks/{directory}/{sample}.bigBed",
    params:
        chrom_sizes=config["genome"]["chromosome_sizes"],
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/bed_to_bigbed/{directory}/{sample}.log",
    shell:
        """
        sort -k1,1 -k2,2n {input.bed} | grep '#' -v | cut -f 1-4 > {input.bed}.tmp &&
        bedToBigBed {input.bed}.tmp {params.chrom_sizes} {output.bigbed} 2> {log} &&
        rm {input.bed}.tmp
        """


localrules:
    save_design,
    make_genomic_bins,
    bed_to_saf,
    validate_peaks,
