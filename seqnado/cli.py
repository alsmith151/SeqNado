import click
import os
import subprocess

FILE = os.path.abspath(__file__)
PACKAGE_DIR = os.path.dirname(FILE)


@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument("method", type=click.Choice(["atac", "chip", "rna", "snp"]))
@click.argument("cookiecutter_options", nargs=-1, type=click.UNPROCESSED)
def cli_config(method, cookiecutter_options, help=False):
    """
    Runs the config for the data processing pipeline.
    """
    cmd = [
        "cookiecutter",
        os.path.join(PACKAGE_DIR, "data/cookiecutter_config", f"config_{method}"),
    ]

    if cookiecutter_options:
        cmd.extend(cookiecutter_options)

    completed = subprocess.run(cmd)


@click.command()
@click.argument("method", type=click.Choice(["atac", "chip", "rna", "snp"]))
@click.argument("files", nargs=-1)
@click.option("-o", "--output", default="design.csv", help="Output file name")
def cli_design(method, files, output="design.csv"):
    """
    Generates a SeqNado design file from a list of files.
    """

    assert len(files) > 0, "No files provided. Please provide a list of files separated by spaces."


    if not method == "chip":
        from seqnado.utils import GenericFastqSamples
        design = GenericFastqSamples.from_files(files).design
    else:
        from seqnado.utils_chipseq import ChipseqFastqSamples
        design = ChipseqFastqSamples.from_files(files).design
    
    design = design.drop(columns=["paired"], errors="ignore")
    design.to_csv(output, index=False)
    

@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument("method", type=click.Choice(["atac", "chip", "rna", "snp"]))
@click.option("-c", "--cores", default=1, help="Number of cores to use", required=True)
@click.option(
    "--preset",
    default="lc",
    help="""Pre-set snakemake job profile to use for pipeline run:
            lc: local conda environment
            ls: local singularity environment
            ss: slurm singularity environment (runs jobs on cluster)
            """,
    type=click.Choice(choices=["lc", "ls", "ss"]),
)
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def cli_pipeline(method, pipeline_options, help=False, cores=1, preset="local"):

    """Runs the data processing pipeline"""

    cmd = [
        "snakemake",
        "-c",
        str(cores),
        "--snakefile",
        os.path.join(PACKAGE_DIR, "workflow", f"snakefile_{method}"),
    ]

    if pipeline_options:
        cmd.extend(pipeline_options)

    if preset == "ss":
        cmd.extend(
            [
                "--profile",
                os.path.abspath(
                    os.path.join(
                        PACKAGE_DIR, "workflow/envs/profiles/profile_drmaa_singularity"
                    )
                ),
            ]
        )
    elif preset == "ls":
        cmd.extend(
            [
                "--profile",
                os.path.abspath(
                    os.path.join(
                        PACKAGE_DIR, "workflow/envs/profiles/profile_singularity"
                    )
                ),
            ]
        )


    with open(f"{PACKAGE_DIR}/data/logo.txt", "r") as f:
        logo = f.read()
    
    print(logo)
    completed = subprocess.run(cmd)
