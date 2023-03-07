import pandas as pd
import itertools
import os
import numpy as np
import re
import pathlib
import subprocess


def get_samplename(path: str):
    p = pathlib.Path(path)
    return re.split(r"_[plus|minus]", p.name)[0]


# Set up details
df = pd.DataFrame(
    itertools.chain.from_iterable([files for files in [snakemake.input.bigwig]]),
    columns=["filename"],
)

df["samplename"] = df["filename"].apply(get_samplename)
df["method"] = df["filename"].apply(lambda x: x.split("/")[-2])
df["strand"] = np.where(df["filename"].str.contains("_plus.bigWig"), "plus", "minus")

file_details = f"{os.path.dirname(snakemake.output.hub)}/hub_details.tsv"
df.set_index("filename").to_csv(file_details, sep="\t")

color_by = snakemake.config["ucsc_hub_details"].get("color_by", None)

if not color_by:
    color_by = ("samplename",)

cmd = " ".join(
    [
        "make-ucsc-hub",
        " ".join(df["filename"]),
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
        " ".join([f"--color-by {c}" for c in color_by]),
        "--group-overlay",
        "samplename",
    ]
)


subprocess.run(cmd, shell=True, check=True)
