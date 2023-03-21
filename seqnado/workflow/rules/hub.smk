import pathlib
import re
import numpy as np


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
        "logs/bed_to_bigbed/{directory}_{sample}.log",
    shell:
        """
        sort -k1,1 -k2,2n {input.bed} > {input.bed}.tmp &&
        bedToBigBed {input.bed}.tmp {params.chrom_sizes} {output.bigbed} || touch {output.bigbed} &&
        rm {input.bed}.tmp
        """

rule generate_hub_for_chipseq_and_atacseq:
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
    container: None,
    script:
        "../scripts/create_hub_genomics.py"


localrules:
    generate_hub_for_chipseq_and_atacseq,
