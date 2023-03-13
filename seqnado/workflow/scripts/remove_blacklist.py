import os
import subprocess

if snakemake.params.blacklist and os.path.exists(snakemake.params.blacklist):

    cmd = f"""
            bedtools intersect -v -b {snakemake.params.blacklist} -a {snakemake.input.bam} > {snakemake.output.bam} &&
            samtools index {snakemake.output.bam} &&
            echo "Removed blacklisted regions" > {snakemake.log}
           """
    subprocess.run(cmd, shell=True, check=True)
else:
    cmd = f"""
           mv {snakemake.input.bam} {snakemake.output.bam} &&
           mv {snakemake.input.bam}.bai {snakemake.output.bam}.bai &&
           echo "No blacklisted regions specified" > {snakemake.log}
           """
    subprocess.run(cmd, shell=True, check=True)
