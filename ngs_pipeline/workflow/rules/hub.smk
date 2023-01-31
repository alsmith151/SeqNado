import pathlib
import re
import numpy as np


rule bed_to_bigbed:
    input:
        bed="peaks/{directory}/{sample}.bed",
    output:
        bigbed="peaks/{directory}/{sample}.bigBed",
    params:
        chrom_sizes=config["genome"]["chromosome_sizes"],
    log:
        "logs/bed_to_bigbed/{directory}_{sample}.log",
    shell:
        """
        sort -k1,1 -k2,2n {input.bed} > {input.bed}.tmp &&
        bedToBigBed {input.bed}.tmp {params.chrom_sizes} {output.bigbed} &&
        rm {input.bed}.tmp
        """


rule generate_hub_for_chipseq_and_atacseq:
    input:
        bigbed=expand(
            "peaks/{method}/{sample}.bigBed",
            method=PEAK_CALL_METHODS
            if ASSAY in ["ChIP", "ATAC"]
            else [
                "",
            ],
            sample=SAMPLE_NAMES_IP if ASSAY == "ChIP" else SAMPLE_NAMES,
        ),
        bigwig=expand(
            "bigwigs/{method}/{sample}.bigWig",
            method=PILEUP_METHODS,
            sample=SAMPLE_NAMES,
        ),
        report="qc/full_qc_report.html",
    output:
        hub=os.path.join(
            config["ucsc_hub_details"]["directory"],
            f"{config['ucsc_hub_details']['name']}.hub.txt",
        ),
    log:
        log=f"logs/{config['ucsc_hub_details']['name']}.hub.log",
    run:
        import pandas as pd
        import itertools

        # Set up details
        df = pd.DataFrame(
            itertools.chain.from_iterable(
                [files for files in [input.bigbed, input.bigwig]]
            ),
            columns=["filename"],
        )

        if ASSAY == "ChIP":
            df[["samplename", "antibody"]] = df["filename"].str.extract(
                r".*/(.*)_(.*)\.(?:bigBed|bigWig)"
            )
        else:
            df["samplename"] = df["filename"].str.extract(r".*/(.*)\.(?:bigBed|bigWig)")

        df["method"] = df["filename"].apply(lambda x: x.split("/")[-2])


        file_details = f"{os.path.dirname(output.hub)}/hub_details.tsv"
        df.set_index("filename").to_csv(file_details, sep="\t")

        color_by = config["ucsc_hub_details"].get("color_by", None)

        if not color_by:
            if df["samplename"].unique().shape[0] == 1:
                color_by = ("antibody",)
            else:
                color_by = ("samplename", "antibody")

        cmd = " ".join(
            [
                "make-ucsc-hub",
                " ".join(df["filename"]),
                "-d",
                file_details,
                "-o",
                os.path.dirname(output.hub),
                "--hub-name",
                config["ucsc_hub_details"]["name"],
                "--hub-email",
                config["ucsc_hub_details"]["email"],
                "--genome-name",
                config["genome"]["name"],
                "--description-html",
                input.report,
                " ".join([f"--color-by {c}" for c in color_by]),
            ]
        )


        shell(cmd)


def get_samplename(path: str):
    p = pathlib.Path(path)
    return re.split(r"_[plus|minus]", p.name)[0]


rule generate_hub_for_rnaseq:
    input:
        bigwig=expand(
            "bigwigs/{method}/{sample}_{strand}.bigWig",
            method=PILEUP_METHODS,
            sample=SAMPLE_NAMES,
            strand=["plus", "minus"],
        ),
        report="qc/full_qc_report.html",
    output:
        hub=os.path.join(
            config["ucsc_hub_details"]["directory"],
            f"{config['ucsc_hub_details']['name']}.hub.txt",
        ),
    log:
        log=f"logs/{config['ucsc_hub_details']['name']}.hub.log",
    run:
        import pandas as pd
        import itertools

        # Set up details
        df = pd.DataFrame(
            itertools.chain.from_iterable([files for files in [input.bigwig]]),
            columns=["filename"],
        )

        df["samplename"] = df["filename"].apply(get_samplename)
        df["method"] = df["filename"].apply(lambda x: x.split("/")[-2])
        df["strand"] = np.where(
            df["filename"].str.contains("_plus.bigWig"), "plus", "minus"
        )

        file_details = f"{os.path.dirname(output.hub)}/hub_details.tsv"
        df.set_index("filename").to_csv(file_details, sep="\t")

        color_by = config["ucsc_hub_details"].get("color_by", None)

        if not color_by:
            color_by = ("samplename",)

        cmd = " ".join(
            [
                "make-ucsc-hub",
                " ".join(df["filename"]),
                "-d",
                file_details,
                "-o",
                os.path.dirname(output.hub),
                "--hub-name",
                config["ucsc_hub_details"]["name"],
                "--hub-email",
                config["ucsc_hub_details"]["email"],
                "--genome-name",
                config["genome"]["name"],
                "--description-html",
                input.report,
                " ".join([f"--color-by {c}" for c in color_by]),
                "--group-overlay",
                "samplename",
            ]
        )


        shell(cmd)


localrules:
    generate_hub_for_chipseq_and_atacseq,
