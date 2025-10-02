import pathlib
import re
import numpy as np


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


rule generate_hub:
    input:
        data=[
            OUTPUTS.bigwig_files,
            OUTPUTS.bigbed_files,
        ],
        report="seqnado_output/seqnado_report.html",
    output:
        hub=OUTPUTS.ucsc_hub.hub_txt,
    log:
        log=f"seqnado_output/logs/{CONFIG.assay_config.ucsc_hub.name}.hub.log",
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    params:
        assay=ASSAY,
        params=CONFIG.assay_config.ucsc_hub,
        has_consensus_peaks=OUTPUTS.has_consensus_peaks,
    script:
        "../scripts/create_hub.py"


localrules:
    generate_hub,
    validate_peaks,
