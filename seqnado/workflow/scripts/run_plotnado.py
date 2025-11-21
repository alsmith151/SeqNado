from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnado.api as pn
import pyranges as pr
import seaborn as sns
from loguru import logger


# Configure logger to write to snakemake log file if available
if "snakemake" in globals():
    logger.remove()
    logger.add(
        snakemake.log[0],
        format="{time} {level} {message}",
        level="DEBUG",
    )
    logger.add(sys.stderr, format="{time} {level} {message}", level="DEBUG")
else:
    logger.add(sys.stderr, format="{time} {level} {message}", level="DEBUG")

logger.info("Starting Plotnado visualization")


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["svg.fonttype"] = "none"

ASSAY = snakemake.params.assay
logger.info(f"Processing {ASSAY} assay")
logger.debug(f"Input files: {snakemake.input.data}")
logger.debug(f"Output plots: {snakemake.output.plots}")
logger.debug(f"Plotting regions: {snakemake.params.regions}")
logger.debug(f"Output directory: {snakemake.params.outdir}")


# Load the tracks into a DataFrame
logger.info("Loading input tracks...")
df = pd.DataFrame([Path(p) for p in snakemake.input.data], columns=["path"])
df["name"] = df["path"].apply(lambda x: x.stem)
df["type"] = df["path"].apply(lambda x: x.suffix)
df["type"] = pd.Categorical(df["type"], categories=[".bigWig", ".bed"], ordered=True)
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
    df["name"] + "-" + df["method"] + "-" + df["normalisation"] + df["type"].astype(str)
)
df["track_name"] = df["track_name"].str.replace("-.", ".")

if ASSAY == "ChIP":
    df["antibody"] = df["name"].str.split("_").str[-1]


# Load the regions
logger.info("Loading plotting regions...")
coords = pr.read_bed(snakemake.params.regions)
logger.info(f"Found {len(coords)} regions to plot")
plotting_format = snakemake.params.plotting_format
logger.info(f"Output format: {plotting_format}")

# Generate the figure
logger.info("Creating Plotnado figure...")
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
        autoscaling_group = None
        style = None
    elif track.type == ".bigWig":
        track_type = "bigwig"
        style = "stairsfilled"

        if ASSAY == "ChIP":
            autoscaling_group = f"{track.antibody}-{track.method}-{track.normalisation}"
        else:
            autoscaling_group = f"{track.method}-{track.normalisation}"

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
        autoscaling_group=autoscaling_group,
    )
    fig.add_track(t)


outdir = Path(snakemake.params.outdir)
logger.info(f"Output directory: {outdir}")
outdir.mkdir(parents=True, exist_ok=True)

for region in coords.df.itertuples():
    fig_name = (
        f"{region.Chromosome}-{region.Start}-{region.End}"
        if not hasattr(region, "Name") and not region.Name
        else region.Name
    )
    region_coords = f"{region.Chromosome}:{region.Start}-{region.End}"
    output_file = outdir / f"{fig_name}.{plotting_format}"
    logger.info(f"Saving plot for region {fig_name}: {output_file}")
    fig.save(output=output_file, gr=region_coords)

logger.info("Saving Plotnado template...")
fig.to_toml(snakemake.output.template)
logger.info("Plotnado visualization complete!")
