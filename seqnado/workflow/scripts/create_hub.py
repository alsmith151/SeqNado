import os
import pandas as pd
import itertools
import numpy as np
import pathlib
import re
from loguru import logger
import tracknado


def get_rna_samplename(path: str):
    p = pathlib.Path(path)
    return re.split(r"_[plus|minus]", p.name)[0]


# Set up details
df = pd.DataFrame(
    snakemake.input.data,
    columns=["fn"],
)


# Use the TrackFiles class to deduplicate files and add metadata
df = tracknado.TrackFiles(files=df, deduplicate=True).files


if snakemake.params.assay == "ChIP":
    df[["samplename", "antibody"]] = df["fn"].str.extract(
        r".*/(.*)_(.*)\.(?:bigBed|bigWig)"
    )
    df["pileup_method"] = df["fn"].apply(lambda x: x.split("/")[-3])
    df['norm'] = df['fn'].apply(lambda x: x.split("/")[-2])

elif snakemake.params.assay == "ATAC":
    df["samplename"] = df["fn"].str.extract(r".*/(.*)\.(?:bigBed|bigWig)")
    df["method"] = df["fn"].apply(lambda x: x.split("/")[-3])
    df["norm"] = df["fn"].apply(lambda x: x.split("/")[-2])

elif snakemake.params.assay == "RNA":
    df["samplename"] = df["fn"].apply(get_rna_samplename)
    df["method"] = df["fn"].apply(lambda x: x.split("/")[-3])
    df["strand"] = np.where(df["fn"].str.contains("_plus.bigWig"), "plus", "minus")
    df["norm"] = df["fn"].apply(lambda x: x.split("/")[-2])

# Create hub design
design = tracknado.TrackDesign.from_design(
    df,
    color_by=snakemake.params.color_by,
    subgroup_by=snakemake.params.subgroup_by
    if any(snakemake.params.subgroup_by)
    else None,
    supergroup_by=snakemake.params.supergroup_by,
    overlay_by=snakemake.params.overlay_by,
)

outdir = pathlib.Path(snakemake.output.hub).parent
hub = tracknado.HubGenerator(
    track_design=design,
    genome=snakemake.params.genome,
    hub_name=snakemake.params.hub_name,
    description_html=pathlib.Path(snakemake.input.report),
    hub_email=snakemake.params.hub_email,
    custom_genome=snakemake.params.custom_genome,
    genome_twobit=snakemake.params.genome_twobit,
    genome_organism=snakemake.params.genome_organism,
    genome_default_position=snakemake.params.genome_default_position,
    outdir=outdir,
)

hub.stage_hub()
design.to_pickle(outdir / ".track_design.pkl")
