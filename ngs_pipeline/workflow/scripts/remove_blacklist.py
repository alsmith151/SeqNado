import os


if snakemake.params.blacklist and os.path.exists(snakemake.params.blacklist):

    cmd = f"""
                        bedtools intersect -v -b {snakemake.params.blacklist} -a {snakemake.input.bam} > {snakemake.input.bam}.tmp &&
                        mv {snakemake.input.bam}.tmp {snakemake.input.bam} &&
                        echo "Removed blacklisted regions" > {snakemake.log}
                        """
    shell(cmd)
else:
    cmd = f"""echo "No blacklisted regions specified" > {log}"""
    shell(cmd)
