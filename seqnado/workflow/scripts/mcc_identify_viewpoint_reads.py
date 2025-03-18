from mcc import mcc
import pathlib



outdir = pathlib.Path(snakemake.params.output_dir)
outdir.mkdir(exist_ok=True, parents=True)
mcc.split_genomic_reads(snakemake.input.bam, snakemake.params.output_dir)