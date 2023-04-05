import os
import pandas as pd
import itertools
import subprocess


# Set up details
df = pd.DataFrame(
    itertools.chain.from_iterable(
        [files for files in [snakemake.input.bigbed, snakemake.input.bigwig]]
    ),
    columns=["filename"],
)

if snakemake.params.assay == "ChIP":
    df[["samplename", "antibody"]] = df["filename"].str.extract(
        r".*/(.*)_(.*)\.(?:bigBed|bigWig)"
    )
else:
    df["samplename"] = df["filename"].str.extract(r".*/(.*)\.(?:bigBed|bigWig)")

df["method"] = df["filename"].apply(lambda x: x.split("/")[-2])


file_details = f"{os.path.dirname(snakemake.output.hub)}/hub_details.tsv"
df.set_index("filename").to_csv(file_details, sep="\t")

color_by = snakemake.config["ucsc_hub_details"].get("color_by", None)

if isinstance(color_by, str):
    color_by = (color_by,)

if not color_by:
    if snakemake.params.assay == "ChIP":
        if df["samplename"].unique().shape[0] == 1:
            color_by = ("antibody",)
        else:
            color_by = ("samplename", "antibody")
    else:
        color_by = ("samplename",)

cmd = [
        "make-ucsc-hub",
        *df["filename"].tolist(),
        "-d",
        file_details,
        "-o",
        os.path.dirname(snakemake.output.hub),
        "--hub-name",
        snakemake.config["ucsc_hub_details"]["name"],
        "--hub-email",
        snakemake.config["ucsc_hub_details"]["email"],
        "--genome-name",
        snakemake.config["genome"]["name"],
        "--description-html",
        snakemake.input.report,
    ]

for cb in color_by:
    cmd.append("--color-by")
    cmd.append(cb)

subprocess.run(cmd)
