import os
import pathlib

import subprocess
import sys

import click
from loguru import logger

from seqnado.helpers import get_genomes

FILE = os.path.abspath(__file__)
PACKAGE_DIR = os.path.dirname(FILE)
GENOMES = get_genomes()


# Init
@click.command(context_settings=dict(ignore_unknown_options=True))
@click.option("--preset", is_flag=True, default=False, help="Use preset genome config")
def cli_init(preset):
    """
    Initializes the seqnado pipeline

    """
    import json
    from loguru import logger
    
    # Get current conda environment
    conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    conda_env_ok = input(
        f"Current conda environment is {conda_env}. Is this correct? (y/n): "
    )
    if conda_env_ok.lower() != "y":
        logger.error(
            "Please activate the correct conda environment and re-run the command"
        )
        sys.exit(1)

    logger.info("Initialising the correct environmental variables for the pipeline")
    os.system(f"bash {PACKAGE_DIR}/init.sh")

    logger.info("Initialising genome files")
    seqnado_config_dir = pathlib.Path("~/.config/seqnado").expanduser()
    seqnado_config_dir.mkdir(parents=True, exist_ok=True)

    genome_template = (
        pathlib.Path(PACKAGE_DIR) / "workflow" / "config" / "genomes_template.json"
    )
    genome_config = seqnado_config_dir / "genome_config.json"
    if not genome_config.exists():
        with open(genome_template, "r") as f:
            genome_template = json.load(f)

        with open(genome_config, "w") as f:
            json.dump(genome_template, f)
    if preset:
        logger.info("Using shared genome config")
        os.system(f"cp {PACKAGE_DIR}/workflow/config/preset_genomes.json" f" {genome_config}")

    logger.info("Initialisation complete")


# Config
@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument("method", type=click.Choice(["atac", "chip", "rna", "snp"]))
@click.option("-r", "--rerun", is_flag=True, help="Re-run the config")
@click.option(
    "-g",
    "--genome",
    default="other",
    help="Genome to use",
    type=str,
)
def cli_config(method, help=False, genome="other", rerun=False):
    """
    Runs the config for the data processing pipeline.
    """
    from importlib.metadata import version
    import seqnado.config as config

    seqnado_version = version("seqnado")

    config.create_config(method, genome, rerun, seqnado_version=seqnado_version)


# Design
@click.command()
@click.argument("method", type=click.Choice(["atac", "chip", "rna", "snp"]))
@click.argument("files", nargs=-1)
@click.option("-o", "--output", default="design.csv", help="Output file name")
def cli_design(method, files, output="design.csv"):
    """
    Generates a SeqNado design file from a list of files.
    """
    import pathlib

    from seqnado.design import Design, DesignIP

    if not files:
        potential_file_locations = [
            ".",
            "fastqs",
            "fastq",
            "data",
            "data/fastqs",
        ]

        for location in potential_file_locations:
            files = list(pathlib.Path(location).glob("*.fastq.gz"))
            if files:
                break

        if not files:
            logger.error("No fastq files provided or found in current directory")
            logger.error(f"""
                         Fastq files can be provided as arguments or found in the following directories:
                         {potential_file_locations}
                         """)
            raise ValueError("No fastq files provided or found in current directory")

    if not method == "chip":
        design = Design.from_fastq_files(files)
    else:
        design = DesignIP.from_fastq_files(files)

    (
        design.to_dataframe()
        .assign(scale_group="all")
        .sort_values("sample_name")
        .to_csv(output, index=False)
    )


# Pipeline
@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument(
    "method",
    type=click.Choice(["atac", "chip", "rna", "snp"]),
)
@click.option("--version", help="Print version and exit", is_flag=True)
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
@click.option(
    "--clean-symlinks",
    is_flag=True,
    help="Remove symlinks created by previous runs. Useful for re-running pipeline after misconfiguration.",
)
@click.option(
    "-s",
    "--scale-resources",
    help="Scale factor the memory and time resources for the pipeline",
    default=1.0,
    type=float,
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase logging verbosity",
)
@click.option(
    "-q",
    "--queue",
    default=None,
    type=str,
    help="Specify the Slurm queue/partition when using the `ss` preset",
)
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def cli_pipeline(
    method,
    pipeline_options,
    help=False,
    preset="local",
    version=False,
    verbose=False,
    clean_symlinks=False,
    scale_resources=1.0,
    queue=None,
):
    """Runs the data processing pipeline"""

    from seqnado.helpers import extract_cores_from_options

    if version:
        from importlib.metadata import version

        _version = version("seqnado")

        print(f"SeqNado version {_version}")
        sys.exit(0)

    if verbose:
        logger.remove()
        logger.add(sys.stderr, level="DEBUG")
    else:
        logger.remove()
        logger.add(sys.stderr, level="INFO")

    pipeline_options, cores = extract_cores_from_options(pipeline_options)

    # Scale the memory and time resources
    os.environ["SCALE_RESOURCES"] = str(scale_resources)

    # Removes old symlinks if requested
    if clean_symlinks:
        logger.info("Cleaning symlinks")
        links = pathlib.Path("seqnado_output/fastqs").glob("*")
        for link in links:
            if link.is_symlink():
                link.unlink()

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
        slurm_profile_path = os.path.abspath(
            os.path.join(
                PACKAGE_DIR, "workflow/envs/profiles/profile_slurm_singularity"
            )
        )
        cmd.extend(["--profile", slurm_profile_path])
        default_resources = [
            f"slurm_partition={queue}" if queue else "slurm_partition=short",
        ]
        cmd.extend(["--default-resources"] + default_resources)

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

    # Print the logo
    with open(f"{PACKAGE_DIR}/data/logo.txt", "r") as f:
        logo = f.read()

    print(logo)

    # Home directory symlinks cause issues with singularity bind mounts
    # to avoid this will change directory to the full resolved path of the current directory
    cwd = str(pathlib.Path(".").resolve())
    os.chdir(cwd)
    os.environ["PWD"] = cwd
    completed = subprocess.run(cmd, cwd=cwd)

    if completed.returncode != 0:
        sys.exit(completed.returncode)
    else:
        sys.exit(0)
