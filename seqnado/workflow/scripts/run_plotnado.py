# %%
import pathlib
import re

import matplotlib.pyplot as plt
import pandas as pd
import plotnado.api as pn
import pyranges as pr
import seaborn as sns
import numpy as np

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["svg.fonttype"] = "none"

# %%
ASSAY = snakemake.params.assay

# %%
# coords = "chr21:28,724,819-29,649,773"

# %%
# bws = list(
#     pathlib.Path(
#         "pytest-current/chipcurrent/2024-11-06_chip_test/seqnado_output/bigwigs/deeptools"
#     ).rglob("*.bigWig")
# )
# bed = [
#     p
#     for p in pathlib.Path(
#         "pytest-current/chipcurrent/2024-11-06_chip_test/seqnado_output/peaks"
#     ).rglob("*.bed")
#     if not "L-tron" in str(p)
# ]

# %%

# Load the tracks into a DataFrame
df = pd.DataFrame(snakemake.input.files, columns=["path"])
df["name"] = df["path"].apply(lambda x: x.stem)
df["type"] = df["path"].apply(lambda x: x.suffix)
df["normalisation"] = np.where(
    df["type"] != ".bed", df["path"].apply(lambda x: x.parts[-2]), ""
)
df["method"] = np.where(
    df["type"] != ".bed",
    df["path"].apply(lambda x: x.parts[-3]),
    df["path"].apply(lambda x: x.parts[-2]),
)
df = df.sort_values(by=["name", "type", "method", "normalisation"])

df["track_name"] = (
    df["name"] + "-" + df["method"] + "-" + df["normalisation"] + df["type"]
)
df["track_name"] = df["track_name"].str.replace("-.", ".")

if ASSAY == "ChIP":
    df["antibody"] = df["name"].str.split("_").str[-1]



# Load the regions
coords = pr.read_bed(snakemake.input.regions)

# Generate the figure
# %%
fig = pn.Figure(
    autospacing=True,
)

fig.add_track(
    "genes",
    genome=snakemake.input.genome,
    gene_style="normal",
    min_gene_length=int(1e3),
    label_y_offset=-75,
    label_loc="right",
    arrow_color="black",
    fontsize=6,
)


names = df["name"].unique()
colors_dict = dict(zip(names, sns.color_palette("tab20", n_colors=len(names))))

for track in df.itertuples():
    if track.type == ".bed":
        track_type = "bed_simple"
        autoscale_group = None
        style = None
    elif track.type == ".bigWig":
        track_type = "bigwig"
        style = "stairsfilled"

        if ASSAY == "ChIP":
            autoscale_group = f"{track.antibody}-{track.method}-{track.normalisation}"
        else:
            autoscale_group = f"{track.method}-{track.normalisation}"

    t = pn.TrackWrapper(
        track_type,
        str(track.path),
        name=track.track_name,
        title=track.track_name,
        color=colors_dict[track.name],
        style=style,
        data_range_style="text",
        data_range_location="right",
        label_on_track=True,
        label_loc="left",
        autoscale_group=autoscale_group,
    )
    fig.add_track(t)

# %%

outdir = pathlib.Path(snakemake.output[0])
for region in coords.df.itertuples():

    fig_name = f"{region.Chromosome}-{region.Start}-{region.End}" if not hasattr(region, "Name") and not region.Name else region.Name
    fig.save(output=outdir / f"{fig_name}.svg", gr=region)
    
# %%
fig.to_toml(snakemake.output[1])




