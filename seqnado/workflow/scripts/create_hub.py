import os
import pandas as pd
import itertools
import numpy as np
import pathlib
import re

from tracknado.design import HubDesign, HubFiles, HubGenerator

def get_rna_samplename(path: str):
    p = pathlib.Path(path)
    return re.split(r"_[plus|minus]", p.name)[0]




# Set up details
df = pd.DataFrame(
    itertools.chain.from_iterable(
        [files for files in [snakemake.input.bigbed, snakemake.input.bigwig]]
    ),
    columns=["fn"],
)

if snakemake.params.assay == "ChIP":
    df[["samplename", "antibody"]] = df["fn"].str.extract(
        r".*/(.*)_(.*)\.(?:bigBed|bigWig)"
    )
    df["method"] = df["fn"].apply(lambda x: x.split("/")[-2])

elif snakemake.params.assay == "ATAC":
    df["samplename"] = df["fn"].str.extract(r".*/(.*)\.(?:bigBed|bigWig)")
    df["method"] = df["fn"].apply(lambda x: x.split("/")[-2])


elif snakemake.params.assay == "RNA":
    df["samplename"] = df["fn"].apply(get_rna_samplename)
    df["method"] = df["fn"].apply(lambda x: x.split("/")[-2])
    df["strand"] = np.where(df["fn"].str.contains("_plus.bigWig"), "plus", "minus")



# Create hub design
design = HubDesign.from_design(df,
                               color_by=snakemake.params.color_by,
                               subgroup_by=snakemake.params.subgroup_by,
                               supergroup_by=snakemake.params.supergroup_by,
                               overlay_by=snakemake.params.overlay_by,
                               ) 



# Create hub

outdir = pathlib.Path(snakemake.output.hub).parent
hub = HubGenerator(track_design=design,
                   genome=snakemake.params.genome,
                   hub_name=snakemake.params.hub_name,
                   hub_shortlabel=snakemake.params.hub_shortlabel,
                   hub_longlabel=snakemake.params.hub_longlabel,
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








                    




