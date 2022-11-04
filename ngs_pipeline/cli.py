import click
import os
import subprocess


@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument("method", type=click.Choice(["atac", "chip", "rna"]))
@click.option("-c", "--cores", default=1, help="Number of cores to use", required=True)
@click.option(
    "--preset",
    default="local-conda",
    help="Pre-set snakemake job profile to use for pipeline run",
    type=click.Choice(choices=["local-conda", "local-singularity","cbrg"])
)
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def cli(method, pipeline_options, help=False, cores=1, preset="local"):

    """Runs the data processing pipeline"""

    file = os.path.abspath(__file__)
    dir_package = os.path.dirname(file)

    if method == "chip":
        cmd = [
            "snakemake",
            "-c",
            str(cores),
            "--snakefile",
            f"{dir_package}/chipseq/snakefile",
        ]
    elif method == "atac":
        cmd = [
            "snakemake",
            "-c",
            str(cores),
            "--snakefile",
            f"{dir_package}/atacseq/snakefile",
        ]
    elif method == "rna":
        cmd = [
            "snakemake",
            "-c",
            str(cores),
            "--snakefile",
            f"{dir_package}/rnaseq/snakefile",
        ]

    if pipeline_options:
        cmd.extend(pipeline_options)
    
    if preset == "cbrg":
        cmd.extend(["--profile", os.path.abspath(os.path.join(dir_package, "profile_drmaa_sigularity"))])

    completed = subprocess.run(cmd)

    if not completed.returncode == 0:
        raise RuntimeError("Pipeline failed. Check the log.")
