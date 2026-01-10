"""Utility functions for SeqNado CLI and general operations."""

from typing import List, Tuple

from loguru import logger

def extract_cores_from_options(options: List[str]) -> Tuple[List[str], int]:
    """
    Extract the number of cores from the snakemake options.
    """

    try:
        cores_flag = options.index("-c")
        cores = int(options[cores_flag + 1])
        options = [
            o for i, o in enumerate(options) if i not in [cores_flag, cores_flag + 1]
        ]
    except ValueError:
        try:
            cores_flag = options.index("--cores")
            cores = int(options[cores_flag + 1])
            options = [
                o
                for i, o in enumerate(options)
                if i not in [cores_flag, cores_flag + 1]
            ]
        except ValueError:
            cores = 1
            logger.warning("No core flag provided. Defaulting to 1 core.")
    except IndexError:
        cores = 1
        options = [o for i, o in enumerate(options) if i not in [cores_flag]]
        logger.warning("Core flag provided but no value given. Defaulting to 1 core.")

    return options, cores


def extract_apptainer_args(options: List[str]) -> Tuple[List[str], str]:
    """
    Extract the apptainer arguments from the snakemake options.
    """
    try:
        apptainer_flag = options.index("--apptainer-args")
        apptainer_args = options[apptainer_flag + 1]
        options = [
            o
            for i, o in enumerate(options)
            if i not in [apptainer_flag, apptainer_flag + 1]
        ]
    except ValueError:
        apptainer_args = ""

    return options, apptainer_args


def remove_unwanted_run_files():
    import glob
    import os
    import shutil

    slurm_files = glob.glob("slurm-*.out")
    sps_files = glob.glob("sps-*")
    simg_files = glob.glob("*.simg")

    for fn in [*slurm_files, *sps_files, *simg_files]:
        try:
            if not os.path.isdir(fn):
                os.remove(fn)
            else:
                shutil.rmtree(fn)

        except Exception as e:
            print(e)


def run_batch_job_on_error(email: str):
    """
    Run a batch job with slurm to send an email on error.
    """

    import subprocess

    slurm_script = f"""#!/bin/bash
    #SBATCH --job-name=seqnado_error_notification
    #SBATCH --mail-type=END
    #SBATCH --mail-user={email}

    echo "An error occurred in the job. Please check the logs for more details."
    """

    with open("error_notification.sh", "w") as f:
        f.write(slurm_script)

    try:
        subprocess.run(["sbatch", "error_notification.sh"], check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to submit slurm job: {e}")


def pepe_silvia():
    print("PEPE SILVIA")
    _pepe_silvia = "https://relix.com/wp-content/uploads/2017/03/tumblr_o16n2kBlpX1ta3qyvo1_1280.jpg"
    return _pepe_silvia