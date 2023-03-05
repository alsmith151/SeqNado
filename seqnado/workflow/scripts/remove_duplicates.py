import os
import subprocess
import seqnado.utils

method = snakemake.config.get("remove_pcr_duplicates_method")

if method == "picard":
    cmd = [
        "picard",
        "MarkDuplicates",
        "-I",
        snakemake.input.bam,
        "-O",
        snakemake.output.bam,
        "-M",
        snakemake.output.metrics,
        "--REMOVE_DUPLICATES true",
        "--CREATE_INDEX true",
        seqnado.utils.check_options(snakemake.params.options),
        ">",
        snakemake.log,
        "2>&1",
    ]
elif method == "samtools":
    cmd = [
        "samtools",
        "rmdup",
        snakemake.input.bam,
        snakemake.output.bam,
        ">",
        snakemake.log,
        "2>&1",
        "&&",
        "samtools",
        "index",
        snakemake.output.bam,
    ]
else:
    cmd = [
        "mv",
        snakemake.input.bam,
        snakemake.output.bam,
        "&&",
        "mv",
        snakemake.input.bam + ".bai",
        snakemake.output.bam + ".bai",
    ]

subprocess.run(" ".join(cmd), shell=True, check=True)


    


