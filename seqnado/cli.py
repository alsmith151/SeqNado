import json
import os
import pathlib
import shutil
import subprocess
import sys

import click
from importlib import resources
from loguru import logger
from seqnado import Assay


FILE = os.path.abspath(__file__)
PACKAGE_DIR = os.path.dirname(FILE)


@click.command(context_settings=dict(ignore_unknown_options=True))
@click.option("--preset", is_flag=True, default=False, help="Use preset genome config")
def cli_init(preset: bool) -> None:
    """
    Initialize SeqNado user environment.

    - Conda is optional (no prompt). If active, we just log it.
    - Uses package resources to read templates and the init script.
    - Creates ~/.config/seqnado/genome_config.json if missing.
    """
    # Conda environments are optional; just log if present
    conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    if conda_env:
        logger.info(f"Conda environment detected: {conda_env}")
    else:
        logger.info("No conda environment detected; proceeding without it.")

    # Initialize apptainer settings if apptainer is available
    if shutil.which("apptainer"):
        try:
            logger.info("Configuring Apptainer/Singularity settings")
            with resources.as_file(resources.files("seqnado").joinpath("init.sh")) as init_path:
                subprocess.run(["bash", str(init_path)], check=True)
        except subprocess.CalledProcessError as e:
            logger.warning(f"Apptainer init script failed (continuing): {e}")
        except Exception as e:
            logger.warning(f"Skipping Apptainer init due to error: {e}")
    else:
        logger.info("Apptainer not found on PATH; skipping container setup.")

    # Prepare genome config path
    seqnado_config_dir = pathlib.Path.home() / ".config" / "seqnado"
    seqnado_config_dir.mkdir(parents=True, exist_ok=True)
    genome_config = seqnado_config_dir / "genome_config.json"

    if genome_config.exists():
        logger.info(f"Found genome config: {genome_config}")
        try:
            genome_data = json.loads(genome_config.read_text())
        except Exception as e:
            logger.warning(f"Could not parse existing genome config ({e}); leaving as-is.")
            genome_data = None

        # Light sanity notice if placeholders are present
        if isinstance(genome_data, dict):
            has_placeholders = any(
                isinstance(v, str) and ("PATH" in v or v == "NA")
                for g in genome_data.values() if isinstance(g, dict)
                for v in g.values()
            )
            if has_placeholders:
                logger.info(
                    "Genome config contains placeholder paths. Update paths before running pipelines."
                )
    else:
        # Write a new genome config from packaged templates
        data_pkg = "seqnado.data"
        template_name = "preset_genomes.json" if preset else "genomes_template.json"
        try:
            with resources.as_file(resources.files(data_pkg).joinpath(template_name)) as tpl_path:
                template = json.loads(pathlib.Path(tpl_path).read_text())
            if preset:
                logger.info("Created genome config from preset genomes.")
            else:
                logger.info(
                    f"Created genome config template at {genome_config}. Please update with valid paths."
                )
            genome_config.write_text(json.dumps(template, indent=4))
        except FileNotFoundError:
            logger.error("Could not locate packaged genome templates; installation may be corrupted.")
            sys.exit(1)
        except Exception as e:
            logger.error(f"Failed to write genome config: {e}")
            sys.exit(1)

    logger.success("Initialization complete")


# Config
@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument("method", type=click.Choice(Assay.all_assay_clean_names()))
@click.option('--dont-make-directories', is_flag=True, help="Do not create directories for the workflow")
@click.option('--all-options', is_flag=True, help="Write all options to the config file, even if they are not used in the workflow")
def cli_config(method, dont_make_directories, all_options):
    """
    Runs the config for the data processing pipeline.
    """
    from importlib.metadata import version
    from seqnado.inputs import Assay
    from seqnado.config.user_input import build_workflow_config, render_config
    from pathlib import Path

    seqnado_version = version("seqnado")
    assay = Assay.from_clean_name(method)
    workflow_config = build_workflow_config(assay, seqnado_version)
    
    if not workflow_config:
        logger.error("Failed to build workflow configuration.")
        sys.exit(1)
    
    if not dont_make_directories:
        dirname = f'{workflow_config.project.date}_{assay.value}_{workflow_config.project.name}'
        outdir = Path(dirname)
        outdir.mkdir(parents=True, exist_ok=True)

        fastq_dir = outdir / "fastqs"
        fastq_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created output directory: {fastq_dir}")

        config_output = outdir / f'config_{assay.value}.yaml'
    
    else:
        config_output = Path(f'config_{assay.value}.yaml')

    # Use packaged template (seqnado.data/config_template.jinja) via importlib.resources
    try:
        with resources.as_file(resources.files("seqnado.data").joinpath("config_template.jinja")) as tpl_path:
            render_config(template=Path(tpl_path), workflow_config=workflow_config, outfile=config_output, all_options=all_options)
    except FileNotFoundError:
        logger.error("Could not locate packaged config template; installation may be corrupted.")
        sys.exit(1)



# Design
@click.command()
@click.argument("method", type=click.Choice(Assay.all_assay_clean_names()))
@click.argument("files", nargs=-1)
@click.option("-o", "--output", default="design.csv", help="Output file name")
@click.option("--merge", is_flag=True, help="Generate a 'merge' column in the design file")
def cli_design(method, files, output="design.csv", merge=False):
    """
    Generates a SeqNado design file from a list of files.
    """
    import pathlib

    from seqnado.inputs import FastqCollection, FastqCollectionForIP

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


    if method in ['chip', 'cat']:
        design = FastqCollectionForIP.from_fastq_files(files)
    else:
        design = FastqCollection.from_fastq_files(files)

    
    df = (
        design.to_dataframe()
                .assign(norm_group="all")
        .sort_values("sample_name")
    )

    if merge and not method == 'rna':
        if 'antibody' in df.columns:
            df['merge'] = df['antibody']
        else:
            df['merge'] = 'consensus'

    df.to_csv(output, index=False)
    logger.success(f"Design file saved to {output}")


# Pipeline
@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument(
    "method",
    required=False,
    type=click.Choice(Assay.all_assay_clean_names()),
)
@click.option("--version", help="Print version and exit", is_flag=True)
@click.option(
    "--preset",
    default="lc",
    help="""Pre-set snakemake job profile to use for pipeline run:
            * lc: local conda environment
            * ls: local singularity environment
            * ss: slurm singularity environment (runs jobs on cluster)
            To use a custom profile, provide the profile options directly using the --profile PROFILE_NAME snakemake argument.
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
    method=None,
    pipeline_options=(),
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

    # If not requesting version, require a method/assay to be provided
    if not method:
        logger.error("No method/assay provided. Provide an assay or use --version to print the version.")
        sys.exit(1)

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
