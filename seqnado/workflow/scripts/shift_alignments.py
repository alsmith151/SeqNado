import subprocess

if snakemake.config.get("shift_atac_reads"):

    cmd = [
        "rsbamtk",
        "shift",
        "-b",
        snakemake.input.bam,
        "-o",
        snakemake.input.bam + ".tmp",
        "&&",
        "samtools",
        "sort",
        snakemake.input.bam + ".tmp",
        "-@",
        str(snakemake.threads),
        "-o",
        snakemake.input.bam,
        "&&",
        "samtools",
        "index",
        snakemake.input.bam,
        "&&",
        "rm",
        snakemake.input.bam + ".tmp",
        "&&",
        f"""echo "Shifted reads" > {snakemake.log}""",
    ]


else:
    cmd = [f'echo "Will not shift reads" > {snakemake.log}']

subprocess.run(" ".join(cmd), shell=True, check=True)
