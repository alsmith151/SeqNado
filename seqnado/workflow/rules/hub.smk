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
            ],
        ),
        "supergroup_by": config["ucsc_hub_details"].get("supergroup_by", None),
    }

    if ASSAY == "RNA":
        hub_params["overlay_by"] = ["samplename", "method"]
        hub_params["subgroup_by"] = ["samplename", "method", "strand"]

    return hub_params


def get_hub_input(wildcards):
    input_files = []

    match ASSAY:
        case "RNA":
            input_files.extend(
                expand(
                    "seqnado_output/bigwigs/{method}/{sample}_{strand}.bigWig",
                    method=config["pileup_method"],
                    sample=SAMPLE_NAMES,
                    strand=["plus", "minus"],
                )
            )

        case "ChIP":
            input_files.extend(
                expand(
                    "seqnado_output/bigwigs/{method}/{sample}.bigWig",
                    method=config["pileup_method"],
                    sample=SAMPLE_NAMES,
                )
            )

            if config["call_peaks"]:
                input_files.extend(
                    expand(
                        "seqnado_output/peaks/{method}/{sample}.bigBed",
                        method=config["peak_calling_method"],
                        sample=SAMPLE_NAMES_IP,
                    )
                )

        case "ATAC":
            input_files.extend(
                expand(
                    "seqnado_output/bigwigs/{method}/{sample}.bigWig",
                    method=config["pileup_method"],
                    sample=SAMPLE_NAMES,
                )
            )

            if config["call_peaks"]:
                input_files.extend(
                    expand(
                        "seqnado_output/peaks/{method}/{sample}.bigBed",
                        method=config["peak_calling_method"],
                        sample=SAMPLE_NAMES,
                    )
                )

        case _:
            input_files = []

    if "merge" in DESIGN.to_dataframe().columns:
        input_files.extend(
            expand(
                "seqnado_output/bigwigs/consensus/{group}.bigWig",
                group=DESIGN.to_dataframe()["merge"].unique(),
            )
        )

    return input_files


def get_peak_files(wildcards):
    peak_files = []

    if config["call_peaks"]:
        if ASSAY == "ChIP":
            peak_files.extend(
                expand(
                    "seqnado_output/peaks/{method}/{sample}.bed",
                    method=config["peak_calling_method"],
                    sample=SAMPLE_NAMES_IP,
                )
            )
        elif ASSAY == "ATAC":
            peak_files.extend(
                expand(
                    "seqnado_output/peaks/{method}/{sample}.bed",
                    method=config["peak_calling_method"],
                    sample=SAMPLE_NAMES,
                )
            )

    if "merge" in DESIGN.to_dataframe().columns:
        peak_files.extend(
            expand(
                "seqnado_output/consensus_peaks/peaks/{group}.bigBed",
                group=DESIGN.to_dataframe().merge.unique(),
            )
        )

    return peak_files


rule save_design:
    output:
        "seqnado_output/design.csv",
    container:
        None
    run:
        DESIGN.to_dataframe().to_csv("seqnado_output/design.csv", index=False)


rule validate_peaks:
    input:
        peaks=get_peak_files,
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
                    if len(peak_entries) < 1:
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
        mem=500,
    log:
        "seqnado_output/logs/bed_to_bigbed/{directory}/{sample}.log",
    shell:
        """
        sort -k1,1 -k2,2n {input.bed} | grep '#' -v > {input.bed}.tmp &&
        bedToBigBed {input.bed}.tmp {params.chrom_sizes} {output.bigbed} &&
        rm {input.bed}.tmp
        """


def get_hub_txt_path():
    import pathlib

    hub_dir = pathlib.Path(config["ucsc_hub_details"]["directory"])
    hub_name = config["ucsc_hub_details"]["name"]
    hub_txt = hub_dir / (f"{hub_name}.hub.txt").replace(" ", "")
    return str(hub_txt)


rule generate_hub:
    input:
        data=get_hub_input,
        report="seqnado_output/qc/alignment_filtered_qc.html",
    output:
        hub=get_hub_txt_path(),
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
