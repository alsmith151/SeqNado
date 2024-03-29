# FAQ

## Pipeline initialisation

### Workflow defines configfile config_chip.yml but it is not present or accessible.

This error occurs when the pipeline is run without a config file present in the working directory. Ensure that seqnado-config has been run before starting the pipeline and that you are in the new directory created by seqnado-config.

Follow the [Pipeline Setup](pipeline.md#create-a-design-file) instructions to create a config file.


## Singularity configuration

### Workflow Error

Failed to pull singularity image from library://asmith151/seqnado/seqnado_pipeline:latest:  
FATAL: Unable to get library client configuration:  
remote has no library client (see https://apptainer.org/docs/user/latest/endpoint.html#no-default-remote)

Fix:

apptainer remote add --no-login SylabsCloud cloud.sylabs.io  
apptainer remote use SylabsCloud  