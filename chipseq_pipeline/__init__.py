import click
import os
import subprocess

file = os.path.abspath(__file__)
dir_package = os.path.dirname(file)


@click.command()
@click.option('-h', '--help', is_flag=True)
@click.argument("mode", type=click.Choice(["make", "show", "clone", "touch"]))
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def cli(mode, pipeline_options, help=False):
    
    '''Runs the data processing pipeline'''   

    cmd = ['python', 
           f'{dir_package}/chipseq_pipeline.py', 
           mode,
           ]
    
    if help:
        cmd.append('--help')
    
    if pipeline_options:
        cmd.extend(pipeline_options)
    
    subprocess.run(cmd)
