from mcc import mcc
import pathlib

outdir = pathlib.Path(snakemake.params.outdir)
outdir.mkdir(exist_ok=True, parents=True)
mcc.identify_ligation_junctions(str(snakemake.input.bam), str(outdir))
