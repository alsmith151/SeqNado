import pyranges as pr
import pandas as pd

chromsizes = (
    pd.read_csv(snakemake.input.chromsizes, sep="\t", header=None).set_index(0)[1].to_dict()
)
genome_tiled = pr.gf.tile_genome(chromsizes, tile_size=snakemake.params.tile_size)
genome_tiled = genome_tiled.df.assign(
    feature="tile", gene_id=lambda df: df.index.astype(str)
).pipe(pr.PyRanges)
genome_tiled.to_gtf(snakemake.output.genome_tiled)
