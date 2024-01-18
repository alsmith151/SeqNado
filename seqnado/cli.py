import click
import os
import subprocess


FILE = os.path.abspath(__file__)
PACKAGE_DIR = os.path.dirname(FILE)


@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument("method", type=click.Choice(["atac", "chip", "rna", "snp"]))
@click.option(
    "-g",
    "--genome",
    default="other",
    help="Genome to use",
    type=click.Choice(
        choices=[
            "dm6",
            "hg19",
            "hg38",
            "hg38_dm6",
            "hg38_mm39",
            "hg38_spikein",
            "mm10",
            "mm39",
            "other",
        ]
    ),
)
def cli_config(method, help=False, genome="other"):
    """
    Runs the config for the data processing pipeline.
    """
    import seqnado.config as config

    config.create_config(method, genome)


@click.command()
@click.argument("method", type=click.Choice(["atac", "chip", "rna", "snp"]))
@click.argument("files", nargs=-1)
@click.option("-o", "--output", default="design.csv", help="Output file name")
def cli_design(method, files, output="design.csv"):
    """
    Generates a SeqNado design file from a list of files.
    """
    import pathlib
    import sys
    from seqnado.utils import Design, DesignIP, FastqFile, FastqFileIP

    if not files:
        files = list(pathlib.Path(".").glob("*.fastq.gz"))

        if not files:
            raise ValueError("No fastq files provided or found in current directory.")

    if not method == "chip":
        design = Design.from_fastq_files([FastqFile(path=fq) for fq in files])
    else:
        from seqnado.utils import DesignIP

        design = DesignIP.from_fastq_files([FastqFileIP(path=fq) for fq in files])

    design.to_dataframe().to_csv(output)


@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument(
    "method",
    type=click.Choice(["atac", "chip", "rna", "snp", "consensus-peaks"]),
)
@click.option("--version", help="Print version and exit", is_flag=True)
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
def cli_pipeline(
    method, pipeline_options, help=False, cores=1, preset="local", version=False
):
    """Runs the data processing pipeline"""

    if version:
        from importlib.metadata import version

        _version = version("seqnado")

        _version = version("seqnado")
        print(f"SeqNado version {_version}")
        return

    cmd = [
        "snakemake",
        "-c",
        str(cores),
        "--snakefile",
        os.path.join(PACKAGE_DIR, "workflow", f"snakefile_{method.replace('-', '_')}"),
    ]

    if pipeline_options:
        cmd.extend(pipeline_options)

    if preset == "ss":
        cmd.extend(
            [
                "--profile",
                os.path.abspath(
                    os.path.join(
                        PACKAGE_DIR, "workflow/envs/profiles/profile_slurm_singularity"
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

    cmd.extend(["--show-failed-logs"])

    with open(f"{PACKAGE_DIR}/data/logo.txt", "r") as f:
        logo = f.read()

    print(logo)

    completed = subprocess.run(cmd)
