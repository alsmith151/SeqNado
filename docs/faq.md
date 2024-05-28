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


## Optional configuration

### Can I merge multiple samples into a single sample?

Yes, you can merge multiple samples into a single sample to generate merged bigWig files and consensus peaks. To do this, you need to create a design file that specifies the samples to be merged. The design file should have a column named "merge" that specifies the samples to be merged e.g.:


| sample | r1 | r2 | deseq2 | merge |
|--------|----|----|--------|-------|
| rna1 | /tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna1_2.fastq.gz | /tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna1_1.fastq.gz | control | control |
| rna2 | /tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna2_2.fastq.gz | /tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna2_1.fastq.gz | control | control |
| rna3 | /tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna3_2.fastq.gz | /tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna3_1.fastq.gz | control | control |
| rna4 | /tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna4_2.fastq.gz | /tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna4_1.fastq.gz | treated | treated |
| rna5 | /tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna5_2.fastq.gz | /tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna5_1.fastq.gz | treated | treated |
| rna6 | /tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna6_2.fastq.gz | /tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna6_1.fastq.gz | treated | treated |
