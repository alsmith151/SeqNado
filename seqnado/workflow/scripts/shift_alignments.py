import subprocess

if snakemake.config.get("shift_atac_reads"):

    cmd = [
        "rsbamtk",
        "shift",
        "-b",
        snakemake.input.bam,
        "-o",
        snakemake.output.bam + ".tmp",
        "&&",
        "samtools",
        "sort",
        snakemake.output.bam + ".tmp",
        "-@",
        str(snakemake.threads),
        "-o",
        snakemake.output.bam,
        "&&",
        f"""echo "Shifted reads" > {snakemake.log}""",
    ]


else:
    cmd = [f'echo "Will not shift reads" > {snakemake.log}', f"mv {snakemake.input.bam} {snakemake.output.bam}"]

subprocess.run(" ".join(cmd), shell=True, check=True)
