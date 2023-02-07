# SeqNado Pipeline

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
    git clone https://github.com/alsmith151/SeqNado.git
    cd SeqNado
    pip install .
    ```

    b) Install from GitHub directly

    ```
    pip install git+https://github.com/alsmith151/SeqNado.git
    ```

1. If you intend to use a cluster e.g. SLURM add the path to the DRMAA interface to your .bashrc:

    ```
    # Access to the DRMAA library: https://en.wikipedia.org/wiki/DRMAA
    echo "export DRMAA_LIBRARY_PATH=/<full-path>/libdrmaa.so" >> ~/.bashrc

    # For CBRG users the command to use is:
    echo "export DRMAA_LIBRARY_PATH=/usr/lib64/libdrmaa.so" >> ~/.bashrc
    ```

## Running the pipeline

1. Setup project directory using seqnado-config

    In the parent directory of the working directory run the following command:

    ```
    seqnado-config atac # ATAC-seq samples
    seqnado-config chip # ChIP-seq/ChIPMentation
    seqnado-config rna # RNA-seq - Not fully tested

    ```

    This will lead you through a series of questions which will create a new project directory, config file and a sample sheet for you to edit.

    cd into the newly made directory and edit the config file and sample sheet.

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
    seqnado atac # ATAC-seq samples
    seqnado chip # ChIP-seq/ChIPMentation
    seqnado rna # RNA-seq - Not fully tested
    ```

    There are several options to visualise which tasks will be performed by the pipeline
    before running.

    ```
    # If using all default settings (this will run on just the login node)
    seqnado atac -c NUMBER_OF_CORES

    # If you want to use the cluster (recommended)
    seqnado atac -c NUMBER_OF_CORES

    # Avoiding network disconnections
    nohup seqnado atac make &
    ```
