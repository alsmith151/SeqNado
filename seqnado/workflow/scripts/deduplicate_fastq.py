import mccnado
import json

stats = mccnado.deduplicate_fastq(snakemake.input.fq1, snakemake.output.fq1, snakemake.input.fq2, snakemake.output.fq2)

with open(snakemake.log[0], "w") as f:
    json.dump(stats, f)