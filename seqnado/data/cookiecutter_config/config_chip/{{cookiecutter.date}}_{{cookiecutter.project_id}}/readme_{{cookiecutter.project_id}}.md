# ChIP-seq Pipeline
Project: {{cookiecutter.project_name}}
Run by: {{cookiecutter.user_name}}
on: {{cookiecutter.date}}

This file contains the instructions for the ChIP-seq pipeline.

## samples

sample sheet: `{{cookiecutter.date}}_{{cookiecutter.project_id}}/samples_{{cookiecutter.project_id}}.csv` can be edited 
 
## config

cookiecutter has generated the file config_chip.yml
The keys marked as essential are required for the pipeline to run.
The keys marked as optional can either be left as the default or adjusted if required.

## run pipeline

The pipeline is run by the following command:

seqnado chip -c N_CORES

To use the singularity container (allows for running the pipeline with a minimal conda environment), 
you will also need to 'bind' paths to the container (this allows for folders outside the current directory to be used i.e. /t1-data).

seqnado chip -c N_CORES --use-singularity --singularity-args "--bind /t1-data --bind /databank "

To run all jobs on the cluster (highly recommended; these options are for slurm i.e. cbrg cluster):

seqnado chip -c N_CORES --drmaa "--cpus-per-task={threads} --mem-per-cpu={resources.mem_mb} --time=24:00:00 "  

Combining both singularity and slurm options:

seqnado chip -c N_CORES --use-singularity --singularity-args "--bind /t1-data --bind /databank " --drmaa "--cpus-per-task={threads} --mem-per-cpu={resources.mem} --time 24:00:00 "

