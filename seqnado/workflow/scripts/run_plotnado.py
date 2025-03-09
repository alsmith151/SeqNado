import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnado.api as pn
import pyranges as pr
import seaborn as sns

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["svg.fonttype"] = "none"

ASSAY = snakemake.params.assay


# Load the tracks into a DataFrame
df = pd.DataFrame([pathlib.Path(p) for p in snakemake.input.data], columns=["path"])
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
coords = pr.read_bed(snakemake.params.regions)
plotting_format = snakemake.params.plotting_format
# Generate the figure
fig = pn.Figure(
    autospacing=True,
)

fig.add_track("scale")

if snakemake.params.genes:
    fig.add_track(
        "genes",
        file=snakemake.params.genes,
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


outdir = pathlib.Path(snakemake.params.outdir)
for region in coords.df.itertuples():
    fig_name = (
        f"{region.Chromosome}-{region.Start}-{region.End}"
        if not hasattr(region, "Name") and not region.Name
        else region.Name
    )
    region_coords = f"{region.Chromosome}:{region.Start}-{region.End}"
    fig.save(output=outdir / f"{fig_name}.{plotting_format}", gr=region_coords)

fig.to_toml(snakemake.output.template)
