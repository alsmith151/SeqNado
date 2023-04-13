# ChIP-seq Pipeline
Project: {{cookiecutter.project_name}}
Run by: {{cookiecutter.user_name}}
on: {{cookiecutter.date}}

This file contains the instructions for the ChIP-seq pipeline.

## samples

sample sheet: `{{cookiecutter.date}}_{{cookiecutter.project_id}}/design.csv` can be edited

## config

cookiecutter has generated the file config_chip.yml
The keys marked as essential are required for the pipeline to run.
The keys marked as optional can either be left as the default or adjusted if required.

## run pipeline

The pipeline is run by the following command:

seqnado chip -c N_CORES --preset 


1. **Copy or link fastq files into the fastq directory**

    Copy:  
    ```cp PATH_TO_FASTQ/example_R1.fastq.gz```

    Symlink: Be sure to use the absolute path for symlinks i.e.  
        ```ln -s /ABSOLUTE_PATH_TO_FASTQ/example_R1.fastq.gz ```  

1. **Set-up sample sheet**

    There are two options for preparing a sample sheet:

    a) Using seqnado-design

    seqnado-design chip fastq/* 


    If samples names match the following conventions then a sample sheet will be generated for your samples:

        ChIP-seq

        * samplename1_Antibody_R1.fastq.gz
        * samplename1_Antibody_R2.fastq.gz
        * samplename1_Input_1.fastq
        * samplename1_Input_2.fastq

    b) Using a custom sample sheet. 

    This is useful for situations in which it can be difficult to appropriately compare IP and Input control samples. 

    * For ChIP-seq samples you will need to create a csv or tsv file with the following columns:

        | sample      | antibody | fq1                              | fq2                              | control              |
        |-------------|----------|----------------------------------|----------------------------------|----------------------|
        | SAMPLE-NAME | ANTIBODY | SAMPLE-NAME_ANTIBODY_R1.fastq.gz | SAMPLE-NAME_ANTIBODY_R2.fastq.gz | CONTROL_SAMPLE_Input |



1. **Running the pipeline**

    All FASTQ files present in the directory will be processed by the pipeline in parallel and
    original FASTQ files will not be modified. If new FASTQ files are added to a pre-run pipeline,
    only the new files will be processed.

    After copying/linking FASTQ files into the working directory and configuring the copy of
    config_chip.yml in the working directory for the current experiment, the pipeline can be run with:

    ```
    seqnado chip # ChIP-seq/ChIPMentation
    ```

    * To visualise which tasks will be performed by the pipeline before running.  
    ```seqnado chip -c 1 --preset ss --dag | dot -Tpng > dag.png```

    * If using all default settings (this will run on just the login node)  
    ```seqnado chip -c NUMBER_OF_CORES```

    * If you want to use the cluster (recommended)  
    ```seqnado chip -c NUMBER_OF_CORES --preset ss```

    * Avoiding network disconnections  
    ```nohup seqnado chip make &```

    **Your processed data can be found in ./seqnado_output**
