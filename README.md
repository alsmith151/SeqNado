
![SeqNado logo](seqnado_logo.jpeg)
# SeqNado Pipeline

Pipeline based on snakemake to process ChIP-seq, ATAC-seq, RNA-seq and short read WGS data for SNP calling.

## Installation

1. Create a basic conda environment (with pip to install python packages) and activate it.  

    ```{bash}
        conda create -n seqnado pip
        conda activate seqnado
    ```

1. Install the pipeline. Three options:
    
    a) Install the package from pip (recommended)

    ```{bash}
        pip install seqnado
    ```

    b) Clone the repositry and install directly.
    ```{bash}
        git clone https://github.com/alsmith151/SeqNado.git
        cd SeqNado
        pip install .
    ```

    c) Install from GitHub directly

    ```{bash}
        pip install git+https://github.com/alsmith151/SeqNado.git
    ```

1. If you intend to use a cluster e.g. SLURM add the path to the DRMAA interface to your .bashrc:

    ```{bash}
        # Access to the DRMAA library: https://en.wikipedia.org/wiki/DRMAA
        echo "export DRMAA_LIBRARY_PATH=/<full-path>/libdrmaa.so" >> ~/.bashrc

        # For CBRG users the command to use is:
        echo "export DRMAA_LIBRARY_PATH=/usr/lib64/libdrmaa.so" >> ~/.bashrc
    ```
  
## Running the pipeline

1. **Setup project directory**

    In the parent directory of desired the working directory run the following command:

    ```
        seqnado-config atac # ATAC-seq samples
        seqnado-config chip # ChIP-seq/ChIPMentation
        seqnado-config rna # RNA-seq - Not fully tested
        seqnado-config snp # snp calling - Not fully tested

    ```

    This will lead you through a series of questions which will create a new project directory, config file and a sample sheet for you to edit.

    cd into the newly made directory and inspect the config file.  

1. **Copy or link fastq files into the fastq directory**

    Copy:  
    ```cp PATH_TO_FASTQ/example_R1.fastq.gz```

    Symlink: Be sure to use the absolute path for symlinks i.e.  
        ```ln -s /ABSOLUTE_PATH_TO_FASTQ/example_R1.fastq.gz ```  

1. **Set-up sample sheet**

    There are two options for preparing a sample sheet:

    a) Using seqnado-design

    ```
        seqnado-design atac fastq/* # ATAC-seq samples
        seqnado-design chip fastq/* # ChIP-seq/ChIPMentation
        seqnado-design rna fastq/* # RNA-seq - Not fully tested
        seqnado-design snp fastq/* # snp calling - Not fully tested

    ```

    If samples names match the following conventions then a sample sheet will be generated for your samples:

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


    b) Using a custom sample sheet. 

    This is useful for situations in which it can be difficult to appropriately compare IP and Input control samples. 

    * For ChIP-seq samples you will need to create a csv or tsv file with the following columns:

        | sample      | antibody | fq1                              | fq2                              | control              |
        |-------------|----------|----------------------------------|----------------------------------|----------------------|
        | SAMPLE-NAME | ANTIBODY | SAMPLE-NAME_ANTIBODY_R1.fastq.gz | SAMPLE-NAME_ANTIBODY_R2.fastq.gz | CONTROL_SAMPLE_Input |



    * For ATAC-seq, RNA-seq or SNP calling samples you will need to create a csv or tsv file with the following columns:

        | sample      | fq1                              | fq2                              |
        |-------------|----------------------------------|----------------------------------|
        | SAMPLE-NAME | SAMPLE-NAME_R1.fastq.gz | SAMPLE-NAME_R2.fastq.gz |


1. **Running the pipeline**

    All FASTQ files present in the directory will be processed by the pipeline in parallel and
    original FASTQ files will not be modified. If new FASTQ files are added to a pre-run pipeline,
    only the new files will be processed.

    After copying/linking FASTQ files into the working directory and configuring the copy of
    config_[*assay*].yml in the working directory for the current experiment, the pipeline can be run with:

    ```
    seqnado atac # ATAC-seq samples
    seqnado chip # ChIP-seq/ChIPMentation
    seqnado rna # RNA-seq - Not fully tested
    seqnado snp # snp calling - Not fully tested
    ```

    * To visualise which tasks will be performed by the pipeline before running.  
    ```seqnado atac -c 1 --preset ss --dag | dot -Tpng > dag.png```

    * If using all default settings (this will run on just the login node)  
    ```seqnado atac -c NUMBER_OF_CORES```

    * If you want to use the cluster (recommended)  
    ```seqnado atac -c NUMBER_OF_CORES --preset ss```

    * Avoiding network disconnections  
    ```nohup seqnado atac make &```

    **Your processed data can be found in ./seqnado_output**
