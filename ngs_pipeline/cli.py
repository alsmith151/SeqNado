import click
import os
import subprocess

FILE = os.path.abspath(__file__)
PACKAGE_DIR = os.path.dirname(FILE)

@click.command()
@click.argument("method", type=click.Choice(["atac", "chip", "rna"]))
def cli_config(method, help=False):
    """
    Runs the config for the data processing pipeline.
    """
    cmd = [
        "cookiecutter",
        os.path.join(PACKAGE_DIR, 'cookiecutter_config', f'config_{method}'),
    ]

    completed = subprocess.run(cmd)

    if not completed.returncode == 0:
        raise RuntimeError("Pipeline config failed. Check the log.")
        
        
@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument("method", type=click.Choice(["atac", "chip", "rna"]))
@click.option("-c", "--cores", default=1, help="Number of cores to use", required=True)
@click.option(
    "--preset",
    default="local-conda",
    help="Pre-set snakemake job profile to use for pipeline run",
    type=click.Choice(choices=["local-conda", "local-singularity", "cbrg"]),
)
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def cli_pipeline(method, pipeline_options, help=False, cores=1, preset="local"):

    """Runs the data processing pipeline"""

    if method == "chip":
        cmd = [
            "snakemake",
            "-c",
            str(cores),
            "--snakefile",
            f"{PACKAGE_DIR}/workflow/snakefile_chip",
        ]
    elif method == "atac":
        cmd = [
            "snakemake",
            "-c",
            str(cores),
            "--snakefile",
            f"{PACKAGE_DIR}/workflow/snakefile_atac",
        ]
    elif method == "rna":
        cmd = [
            "snakemake",
            "-c",
            str(cores),
            "--snakefile",
            f"{PACKAGE_DIR}/workflow/snakefile_rna",
        ]

    if pipeline_options:
        cmd.extend(pipeline_options)

    if preset == "cbrg":
        cmd.extend(
            [
                "--profile",
                os.path.abspath(os.path.join(PACKAGE_DIR, "profile_drmaa_sigularity")),
            ]
        )

    completed = subprocess.run(cmd)

    if not completed.returncode == 0:
        raise RuntimeError("Pipeline failed. Check the log.")
