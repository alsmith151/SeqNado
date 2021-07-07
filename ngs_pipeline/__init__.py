import click
import os
import subprocess

file = os.path.abspath(__file__)
dir_package = os.path.dirname(file)


@click.command(context_settings=dict(ignore_unknown_options=True))
@click.option('-h', '--help', is_flag=True)
@click.argument("method", type=click.Choice(['atac', 'chip', 'rna']))
@click.argument("mode", type=click.Choice(["make", "show", "clone", "touch"]))
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def cli(method, mode, pipeline_options, help=False):
    
    '''Runs the data processing pipeline'''   

    if method == 'chip' or method == 'atac':
        cmd = ['python', f'{dir_package}/pipeline_atac_chipseq.py', mode,]
    else:
        cmd = ['python', f'{dir_package}/pipeline_rna.py', mode]


    
    if help:
        cmd.append('--help')
    
    if pipeline_options:
        cmd.extend(pipeline_options)
    
    subprocess.run(cmd)
