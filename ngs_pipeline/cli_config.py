import click
import os
import subprocess


@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument("method", type=click.Choice(["atac", "chip", "rna"]))
def cli_config(method, help=False):
    """
    Runs the config for the data processing pipeline.
    """
    file = os.path.abspath(__file__)
    dir_package = os.path.dirname(file)

    cmd = [
        "cookiecutter",
        os.path.join(dir_package, 'cookiecutter_config', f'config_{method}'),
    ]

    completed = subprocess.run(cmd)

    if not completed.returncode == 0:
        raise RuntimeError("Pipeline config failed. Check the log.")
