from mcc.cooler import CoolerMerger, CoolerBinsLinker


infiles = snakemake.input.mcools
outfile = snakemake.output.mcool


# Merge mcool files
CoolerMerger(infiles, outfile).merge_files()

# Link bins
CoolerBinsLinker(outfile).link_bins()