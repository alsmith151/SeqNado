# Pipeline

## Configuration

The pipeline is configured using a YAML file: e.g. `config_atac.yml`, `config_chip.yml`. We highly recommend using the `seqnado-config` command to generate the configuration file as this will prompt the user for the required information and ensure that the configuration file is valid. The configuration file can be edited manually if required e.g. using `nano` or the VS Code text editor.

### Generate the working directory and configuration file

The following command will generate the working directory and configuration file for the ATAC-seq pipeline:

```bash
seqnado-config atac
```

You should get somthing like this:

```bash
$ seqnado-config chip
  [1/23] user_name (Your name): asmith
  [2/23] Select date
    1 - 2024-01-13
    Choose from [1] (1):
  [3/23] project_name (Project name): TEST
  [4/23] Select project_id
    1 - test
    Choose from [1] (1): 1
  [5/23] genome (hg38):
  [6/23] chromosome_sizes (/ceph/project/milne_group/shared/seqnado_reference/hg38/UCSC/sequence/hg38.chrom.sizes):
  [7/23] indicies (/ceph/project/milne_group/shared/seqnado_reference/hg38/UCSC/bt2_index/hg38):
  [8/23] gtf (/ceph/project/milne_group/shared/seqnado_reference/hg38/UCSC/genes/hg38.ncbiRefSeq.gtf):
  [12/23] Select remove_pcr_duplicates_method
    1 - picard
    2 - deeptools
    Choose from [1/2] (1): 1
  [13/23] Select remove_blacklist
    1 - yes
    2 - no
    Choose from [1/2] (1): 1
  [14/23] blacklist (/ceph/project/milne_group/shared/seqnado_reference/hg38/hg38-blacklist.v2.bed.gz):
  [15/23] Select make_bigwigs
    1 - yes
    2 - no
    Choose from [1/2] (1): 1
  [16/23] Select pileup_method
    1 - deeptools
    2 - homer
    Choose from [1/2] (1): 1
  [17/23] Select make_heatmaps
    1 - yes
    2 - no
    Choose from [1/2] (1): 1
  [18/23] Select call_peaks
    1 - yes
    2 - no
    Choose from [1/2] (1): 1
  [19/23] Select peak_calling_method
    1 - macs
    2 - lanceotron
    3 - homer
    Choose from [1/2/3] (1): 2
  [20/23] Select make_ucsc_hub
    1 - yes
    2 - no
    Choose from [1/2] (1): 1
  [21/23] UCSC_hub_directory (path/to/ publically accessible location on the server): /project/milne_group/datashare/asmith/chipseq/TEST_HUB
  [22/23] email (Email address (UCSC required)): alastair.smith@ndcls.ox.ac.uk
  [23/23] Select color_by
    1 - samplename
    2 - method
    Choose from [1/2] (1): 1
```

This will generate the following files:

```bash
$ tree 2024-01-13_test/

2024-01-13_test/
├── config_chip.yml
└── readme_test.md

0 directories, 2 files
```

### Edit the configuration file (if required)

The configuration file can be edited manually if required e.g. using `nano` or the VS Code text editor. Use this if you have made an error in the configuration file or if you want to change it for any other reason.

!!! Warning
    If you edit the configuration file manually, you must ensure that it is valid YAML syntax (ensure that you do not delete any colons, commas, or change the indentation). You can check that the file is valid using the following command:

```bash
nano config_chip.yml # Note to exit nano press ctrl+x and then "y" followed by "enter" to save
```

### Create a design file (optional)

#### Infer sample names from fastq file names

If the fastq files are named in a way that seqnado can infer the sample names, then a design file will be generated automatically:

  ChIP-seq

  * samplename1_Antibody_R1.fastq.gz
  * samplename1_Antibody_R2.fastq.gz
  * samplename1_Input_1.fastq.gz
  * samplename1_Input_2.fastq.gz

  For ATAC-seq:

  * sample-name-1_R1.fastq.gz
  * sample-name-1_R2.fastq.gz
  * sample-name-1_1.fastq.gz
  * sample-name-1_2.fastq.gz

  For RNA-seq:

  * sample-name-1_R1.fastq.gz
  * sample-name-1_R2.fastq.gz
  * sample-name-1_1.fastq.gz
  * sample-name-1_2.fastq.gz


#### Use `seqnado-design` to generate a design file

If the fastq files are not named in a way that seqnado can infer the sample names, then a design file can be generated using the `seqnado-design` command. You'll need to enter the working directory and generate a design file:

```bash
cd 2024-01-13_test/
seqnado-design chip /path/to/fastq/files/* # Note that you can use tab completion to complete the path to the fastq files
```

This will generate a design file called `design.csv` in the working directory.


#### ATAC|RNA-seq design file

An ATAC-seq or RNA-seq design file should look something like this:

```bash
,r1,r2
rna,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna_2.fastq.gz,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna_1.fastq.gz
```

!!! Note
    The design file is a CSV file with the following columns:
      * The first column is the sample name
      * `r1` - The path to the read 1 fastq file
      * `r2` - The path to the read 2 fastq file


#### ChIP-seq design file

A ChIP assay design file should look something like this:

```bash
,ip_r1,ip_r2,control_r1,control_r2,ip,control
CTCF,CTCF_CTCF_2.fastq.gz,CTCF_CTCF_1.fastq.gz,CTCF_input_2.fastq.gz,CTCF_input_1.fastq.gz,CTCF,input
```

!!! Note
    The design file is a CSV file with the following columns:
      * The first column is the sample name
      * `ip_r1` - The path to the IP read 1 fastq file
      * `ip_r2` - The path to the IP read 2 fastq file
      * `control_r1` - The path to the control read 1 fastq file
      * `control_r2` - The path to the control read 2 fastq file
      * `ip` - The name of the IP sample
      * `control` - The name of the control sample


#### RNA-seq design file

An RNA-seq design file should look something like this:

```bash
,r1,r2
rna,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna_2.fastq.gz,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna_1.fastq.gz
```

If you want to run DeSeq2, then you will need to add an additional column to the design file to indicate which samples are in the control group:

```bash
,r1,r2,deseq2
rna1,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna1_2.fastq.gz,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna1_1.fastq.gz,control
rna2,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna2_2.fastq.gz,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna2_1.fastq.gz,control
rna3,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna3_2.fastq.gz,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna3_1.fastq.gz,control
rna4,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna4_2.fastq.gz,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna4_1.fastq.gz,treated
rna5,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna5_2.fastq.gz,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna5_1.fastq.gz,treated
rna6,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna6_2.fastq.gz,/tmp/pytest-of-asmith/pytest-7/data2/2024-01-13_rna_test/rna6_1.fastq.gz,treated
```


## Running the pipeline

### Ensure files are in the correct location

Before running the pipeline, ensure that the fastq files, and design are in the correct location:

```bash
# Fastq files
ln -s /path/to/fastq_files/ /path/to/working-directory/made-by-seqnado-config/

# Design
mv /path/to/design.csv /path/to/working-directory/made-by-seqnado-config/
```

### Check you are in the correct directory

```bash
$ ls -l
-rw-r--r-- 1 asmith asmithgrp    1845 Jan 13 10:50 config_rna.yml
-rw-r--r-- 1 asmith asmithgrp   14784 Jan 13 10:50 deseq2_test.qmd
-rw-r--r-- 1 asmith asmithgrp     155 Jan 13 14:40 design.csv
-rw-r--r-- 1 asmith asmithgrp 3813176 Jan 13 10:50 rna_1.fastq.gz
-rw-r--r-- 1 asmith asmithgrp 3836966 Jan 13 10:50 rna_2.fastq.gz
```


### Ensure that the pipeline will not stop when you log out

```bash
tmux new -s NAME_OF_SESSION
```


### Check you have activated the conda environment

```bash
conda activate seqnado
```

### Run the pipeline

The pipeline can be run using the following command:

```bash
seqnado [atac|chip|rna|snp] -c <number of cores> --preset [ss|ls]
```

An actual example would be:

```bash
seqnado rna -c 8 --preset ss
```

!!! Note
    * To visualise which tasks will be performed by the pipeline before running.
    ```seqnado atac -c 1 --preset ss --dag | dot -Tpng > dag.png```


## Pipeline Errors

Check the log file for errors:

```bash
# Look at all of the log files
cat seqnado_error.log

# Look for errors in the log files
cat seqnado_error.log | grep exception -A 10 -B 10
```
