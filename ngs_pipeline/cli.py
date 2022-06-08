import click
import os
import subprocess

@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument("method", type=click.Choice(['atac', 'chip', 'rna']))
@click.option("-c", "--cores", default=1, help="Number of cores to use")
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def cli(method, pipeline_options, help=False, cores=1):
    
    '''Runs the data processing pipeline'''   

    file = os.path.abspath(__file__)
    dir_package = os.path.dirname(file)

    if method == 'chip':
        cmd = ['snakemake', "-c", str(cores), "--snakefile", f'{dir_package}/chipseq/snakefile',]
    else:
        cmd = ['python', f'{dir_package}/pipeline_rna.py']
    
    if pipeline_options:
        cmd.extend(pipeline_options)
    
    completed = subprocess.run(cmd)

    if not completed.returncode == 0:
        raise RuntimeError("Pipeline failed. Check the log.")
