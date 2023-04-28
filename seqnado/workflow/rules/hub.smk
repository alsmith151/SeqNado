import pathlib
import re
import numpy as np

def get_hub_params(config):
    hub_params = {
        "hub_name": config["ucsc_hub_details"]["name"],
        "hub_short_label": config["ucsc_hub_details"].get("short_label", config["ucsc_hub_details"]["name"]),
        "hub_long_label": config["ucsc_hub_details"].get("long_label", config["ucsc_hub_details"]["name"]),
        "hub_email": config["ucsc_hub_details"]["email"],
        "genome": config["genome"]["name"],
        "custom_genome": config["genome"].get("custom_genome", False),
        "genome_twobit": config["genome"].get("twobit"),
        "genome_organism": config["genome"].get("organism"),
        "genome_default_pos": config["genome"].get("default_pos", "chr1:1-1000000"),
        "color_by": config["ucsc_hub_details"].get("color_by", ["samplename", ]),
        "overlay_by": config["ucsc_hub_details"].get("overlay_by", None),
        "subgroup_by": config["ucsc_hub_details"].get("subgroup_by", ["method",]),
        "supergroup_by": config["ucsc_hub_details"].get("supergroup_by", None),
    }

    if ASSAY == "RNA":
        hub_params["overlay_by"] = ["samplename", "strand"]
        hub_params["subgroup_by"] = ["method", "strand"]
    

    return hub_params



rule bed_to_bigbed:
    input:
        bed="seqnado_output/peaks/{directory}/{sample}.bed",
    output:
        bigbed="seqnado_output/peaks/{directory}/{sample}.bigBed",
    params:
        chrom_sizes=config["genome"]["chromosome_sizes"],
    resources:
        mem_mb=500
    log:
        "seqnado_output/logs/bed_to_bigbed/{directory}_{sample}.log",
    shell:
        """
        sort -k1,1 -k2,2n {input.bed} | grep '#' -v > {input.bed}.tmp &&
        bedToBigBed {input.bed}.tmp {params.chrom_sizes} {output.bigbed} &&
        rm {input.bed}.tmp
        """

rule generate_hub:
    input:
        bigbed=expand(
            "seqnado_output/peaks/{method}/{sample}.bigBed",
            method=PEAK_CALL_METHODS
            if ASSAY in ["ChIP", "ATAC"]
            else [
                "",
            ],
            sample=SAMPLE_NAMES_IP if ASSAY == "ChIP" else SAMPLE_NAMES,
        ),
        bigwig=expand(
            "seqnado_output/bigwigs/{method}/{sample}.bigWig",
            method=PILEUP_METHODS,
            sample=SAMPLE_NAMES,
        ),
        report="seqnado_output/qc/full_qc_report.html",
    output:
        hub=os.path.join(
            config["ucsc_hub_details"]["directory"],
            f"{config['ucsc_hub_details']['name']}.hub.txt",
        ),
    log:
        log=f"seqnado_output/logs/{config['ucsc_hub_details']['name']}.hub.log",
    params:
        assay = ASSAY,
        **get_hub_params(config),
    script:
        "../scripts/create_hub.py"


localrules:
    generate_hub
