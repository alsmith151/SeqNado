import pathlib
import re
import numpy as np



def get_hub_params(config):
    hub_params = {
        "hub_name": config["ucsc_hub_details"]["name"],
        "hub_email": config["ucsc_hub_details"]["email"],
        "genome": config["genome"]["name"],
        "custom_genome": config["genome"].get("custom_genome", False),
        "genome_twobit": config["genome"].get("twobit"),
        "genome_organism": config["genome"].get("organism"),
        "genome_default_position": config["genome"].get(
            "default_pos", "chr1:1-1000000"
        ),
        "color_by": config["ucsc_hub_details"].get(
            "color_by",
            [
                "samplename",
            ],
        ),
        "overlay_by": config["ucsc_hub_details"].get("overlay_by", None),
        "subgroup_by": config["ucsc_hub_details"].get(
            "subgroup_by",
            [
                "method",
                "norm"
            ],
        ),
        "supergroup_by": config["ucsc_hub_details"].get("supergroup_by", "ext"),
    }

    if ASSAY == "RNA":
        hub_params["overlay_by"] = ["samplename", "method", "norm"]
        hub_params["subgroup_by"] = ["method", "norm", "strand"]

    return hub_params


rule save_design:
    output:
        "seqnado_output/design.csv",
    container:
        None
    run:
        DESIGN.to_dataframe().to_csv("seqnado_output/design.csv", index=False)


rule validate_peaks:
    input:
        peaks=OUTPUT.peaks,
    output:
        sentinel="seqnado_output/peaks/.validated",
    container:
        None
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
        mem="1GB",
    log:
        "seqnado_output/logs/bed_to_bigbed/{directory}/{sample}.log",
    shell:
        """
        sort -k1,1 -k2,2n {input.bed} | grep '#' -v > {input.bed}.tmp &&
        bedToBigBed {input.bed}.tmp {params.chrom_sizes} {output.bigbed} 2> {log} &&
        rm {input.bed}.tmp
        """


rule generate_hub:
    input:
        data=[
            OUTPUT.bigwigs,
            OUTPUT.bigbed,
        ],
        report="seqnado_output/qc/alignment_filtered_qc.html",
    output:
        hub=OUTPUT.ucsc_hub.hub_txt,
    log:
        log=f"seqnado_output/logs/{config['ucsc_hub_details']['name']}.hub.log".strip(),
    container:
        None
    params:
        **get_hub_params(config),
        assay=ASSAY,
        has_consensus_peaks="merge" in DESIGN.to_dataframe().columns,
    script:
        "../scripts/create_hub.py"


localrules:
    generate_hub,
    validate_peaks,
