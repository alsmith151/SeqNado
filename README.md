# NGS Pipeline

Pipeline based on snakemake to process ChIP-seq, ATAC-seq and RNA-seq data.

## Installation

1. Create a basic conda environment (with pip to install python packages) and activate it.

```{bash}
conda create -n np pip
conda activate np
````

1. Install the pipeline. Two options:

a) Clone the repositry and install directly

```
git clone https://github.com/alsmith151/ngs_pipeline.git
cd ngs_pipeline
pip install .
```

b) Install from GitHub directly

```
pip install git+https://github.com/alsmith151/ngs_pipeline.git
```

1. If you intend to use a cluster e.g. SLURM add the path to the DRMAA interface to your .bashrc:

```
# Access to the DRMAA library: https://en.wikipedia.org/wiki/DRMAA
echo "export DRMAA_LIBRARY_PATH=/<full-path>/libdrmaa.so" >> ~/.bashrc

# For CBRG users the command to use is:
echo "export DRMAA_LIBRARY_PATH=/usr/lib64/libdrmaa.so" >> ~/.bashrc
```

## Running the pipeline

1. Create a working directory


To run the pipeline you will need to create a working directory for the pipeline run:

```
mkdir RS411_EPZ5676/
cd RS411_EPZ5676/
```

The pipeline will be executed here and all files will be generated
in this directory.

1. Get and edit the pipeline configuration file.

The configuration file [config_X.yml](https://github.com/alsmith151/ngs_pipeline/blob/master/config_atac.yml) enables 
parameterisation of the pipeline run with user specific settings. Furthermore,
it also provides paths to essential files for the pipeline run (e.g., bowtie2 indices).
The paths supplied do not have to be in the same directory as the pipeline.

A copy of config_X.yml can be downloaded from GitHub using:
```
wget https://raw.githubusercontent.com/alsmith151/ngs_pipeline/master/config_atac.yml # ATAC
wget https://raw.githubusercontent.com/alsmith151/ngs_pipeline/master/config_chip.yml # ChIP
wget https://raw.githubusercontent.com/alsmith151/ngs_pipeline/master/config_rna.yml  # RNA

```

This [yaml](https://yaml.org/spec/1.2/spec.html) file can be edited using standard text editors e.g.

```
nano config.yml
```

1. Copy or link fastq files into the working directory

Copy:

```
cp PATH_TO_FASTQ/example_R1.fastq.gz
```

Symlink:

```
# Be sure to use the absolute path for symlinks
ln -s /ABSOLUTE_PATH_TO_FASTQ/example_R1.fastq.gz
```

1. Set-up sample sheet.

There are two options for preparing a sample sheet:

a) Using sample naming. If samples names match the following conventions then a sample sheet will be generated for your samples:

ChIP-seq

* samplename1_Antibody_R1.fastq.gz
* samplename1_Antibody_R2.fastq.gz
* samplename1_Input_1.fastq
* samplename1_Input_2.fastq

For ATAC-seq:

* sample-name-1_R1.fastq.gz
* sample-name-1_R2.fastq.gz
* sample-name-1_1.fastq
* sample-name-1_2.fastq

For RNA-seq:

* sample-name-1_R1.fastq.gz
* sample-name-1_R2.fastq.gz
* sample-name-1_1.fastq
* sample-name-1_2.fastq

b) Using a custom sample sheet. This is useful for situations in which it can be difficult to appropriately compare IP and Input control samples. For ChIP-seq samples you will need to create a csv or tsv file with the following columns:

| sample      | antibody | fq1                              | fq2                              | control              |
|-------------|----------|----------------------------------|----------------------------------|----------------------|
| SAMPLE-NAME | ANTIBODY | SAMPLE-NAME_ANTIBODY_R1.fastq.gz | SAMPLE-NAME_ANTIBODY_R2.fastq.gz | CONTROL_SAMPLE_Input |


1. Running the pipeline

All FASTQ files present in the directory will be processed by the pipeline in parallel and
original FASTQ files will not be modified. If new FASTQ files are added to a pre-run pipeline,
only the new files will be processed.

After copying/linking FASTQ files into the working directory and configuring the copy of
config.yml in the working directory for the current experiment, the pipeline can be run with:

```
ngs-pipeline atac # ATAC-seq samples
ngs-pipeline chip # ChIP-seq/ChIPMentation
ngs-pipeline rna # RNA-seq - Not fully tested 
```

There are several options to visualise which tasks will be performed by the pipeline
before running. 

```
# If using all default settings and using a cluster
ngs-pipeline atac -c NUMBER_OF_CORES

# If not using a cluster, run in local mode.
ngs-pipeline atac 

# Avoiding network disconnections
nohup ngs-pipeline atac make &
```








