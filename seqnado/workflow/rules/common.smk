
rule save_design:
    output:
        "seqnado_output/design.csv",
    container:
        None
    run:
        DESIGN.to_dataframe().to_csv("seqnado_output/design.csv", index=False)


rule make_genomic_bins:
    input:
        chrom_sizes=config["genome"]["chromosome_sizes"],
    params:
        bin_size=config["bin_size"],
    output:
        bed="seqnado_output/resources/genomic_bins.bed",
    log:
        "seqnado_output/logs/genomic_bins.log",
    container:
        None
    shell:
        """
        cooler makebins {input.chrom_sizes} {params.bin_size} > {output.bed}
        """
