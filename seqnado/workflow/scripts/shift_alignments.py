import subprocess

shift_reads = snakemake.config.get("shift_atac_reads", False)

if shift_reads:

    cmd_shift = " ".join(["rsbamtk", "shift", "-b", snakemake.input.bam, "-o", snakemake.output.bam + ".tmp"])
    cmd_sort = " ".join(["samtools", "sort", snakemake.output.bam + ".tmp", "-@", str(snakemake.threads), "-o", snakemake.output.bam])
    cmd_index = " ".join(["samtools", "index", snakemake.output.bam])
    cmd_log = " ".join([f'echo "Shifted reads" > {snakemake.log}'])

    cmd = " && ".join([cmd_shift, cmd_sort, cmd_index, cmd_log])


else:
    cmd = [f'echo "Will not shift reads" > {snakemake.log}', 
           f"mv {snakemake.input.bam} {snakemake.output.bam}",
           f"mv {snakemake.input.bam}.bai {snakemake.output.bam}.bai"]
    
    cmd = " && ".join(cmd)

subprocess.run(cmd, shell=True, check=True)
