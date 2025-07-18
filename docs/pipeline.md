# Pipeline

## Configuration

The pipeline is configured using a YAML file: e.g. `config_atac.yml`, `config_chip.yml`. We highly recommend using the `seqnado-config` command to generate the configuration file as this will prompt the user for the required information and ensure that the configuration file is valid. The configuration file can be edited manually if required e.g. using `nano` or the VS Code text editor.

### Generate the working directory and configuration file

The following command will generate the working directory and configuration file for the ATAC-seq pipeline:

```bash
seqnado-config chip

# options 
-r, --rerun # Re-runs the config in existing seqnado directory

```

You should get something like this:

```bash
$ seqnado-config chip
  What is your project name? [cchahrou_project]: cchahrou_project
  What is the genome? [hg38]: hg38
  Perform fastqscreen? (yes/no) [no]: yes
  Path to fastqscreen config: [/PATH/TO/fastq_screen.conf]: /PATH/TO/fastq_screen.conf
  Do you want to remove blacklist regions? (yes/no) [yes]: yes
  Remove PCR duplicates? (yes/no) [yes]: yes
  Remove PCR duplicates method: [picard/samtools]: picard
  Calculate library complexity? (yes/no) [no]: yes
  Do you have spikein? (yes/no) [no]: yes
  Normalisation method: [orlando/with_input]: orlando
  Reference genome: [hg38]: hg38
  Spikein genome: [dm6]: dm6
  Do you want to make bigwigs? (yes/no) [no]: yes
  Pileup method: [deeptools/homer]: deeptools
  Do you want to make heatmaps? (yes/no) [no]: yes
  Do you want to call peaks? (yes/no) [no]: yes
  Generate consensus counts from Design merge column? (yes/no) [no]: yes
  Peak caller: [lanceotron/macs/homer/seacr]: lanceotron
  Do you want to make a UCSC hub? (yes/no) [no]: yes
  UCSC hub directory: [seqnado_output/hub/]: seqnado_output/hub/
  What is your email address? [email@example.com]: email for UCSC
  Color by (for UCSC hub): [samplename]: samplename
  Generate GEO submission files (MD5Sums, read count summaries...)? (yes/no) [no]: yes
  Perform plotting? (yes/no) [no]: yes
  Path to bed file with coordinates for plotting [None]: 
  Path to bed file with genes. [None]: 
  Directory '2025-02-11_chip_cchahrou_project/' has been created with the 'config_chip.yml' file.
```

This will generate the following files:

```bash
$ tree 2025-02-11_chip_cchahrou_project/

2025-02-11_chip_cchahrou_project/
├── config_chip.yml
└── fastq/

1 directory, 1 file
```

### Check the tool options

In the newly created config yaml file, check the tool options **especially for rna quantification**!

### Edit the configuration file (if required)

The configuration file can be edited manually if required e.g. using `nano` or the VS Code text editor. Use this if you have made an error in the configuration file or if you want to change it for any other reason.

!!! Warning
    If you edit the configuration file manually, you must ensure that it is valid YAML syntax (ensure that you do not delete any colons, commas, or change the indentation). 

  ```bash
  nano config_chip.yml # Note to exit nano press ctrl+x and then "y" followed by "enter" to save
  ```

## Organise your fastq files

use symlinks from your raw_data:

```bash
cd 2025-02-11_chip_cchahrou_project/fastq

ln -s /path/to/raw_data/SampleName_S1_L001_R1_001.fastq.gz samplename1_Antibody_R1.fastq.gz
ln -s /path/to/raw_data/SampleName_S1_L001_R2_001.fastq.gz samplename1_Antibody_R2.fastq.gz
```

### Infer sample names from fastq file names

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


## Design file

If the fastq files are not named in a way that seqnado can infer the sample names, then a design file can be generated using the `seqnado-design` command. You'll need to enter the working directory and generate a design file:


```bash
cd ..
seqnado-design chip fastq/* 
# Note that you can use tab completion to complete the path to the fastq files
```

This will generate a design file called `design.csv` in the working directory.

!!! Warning
    You need to specify the fastq files in the command line to use for the design generation e.g. in the current working directory:  
    ```bash 
    seqnado-design chip *.fastq.gz
    ```

### Merging replicates or samples

To merge samples for counting or bigwig/peak generation add a consensus_group column to the design file

```bash
sample_name,r1,r2,norm_group,consensus_group
atac,/ceph/project/milne_group/cchahrou/software/SeqNado/2025-02-11_chip_cchahrou_project/atac_1.fastq.gz,/ceph/project/milne_group/cchahrou/software/SeqNado/2025-02-11_chip_cchahrou_project/atac_2.fastq.gz,all,consensus_group
atac2,/ceph/project/milne_group/cchahrou/software/SeqNado/2025-02-11_chip_cchahrou_project/atac_1.fastq.gz,/ceph/project/milne_group/cchahrou/software/SeqNado/2025-02-11_chip_cchahrou_project/atac_2.fastq.gz,all,consensus_group
```

This will merge both to make a `consensus_group` bigwig and peak file

Consensus counts can be made from the merged peaks when consensus_counts is True in config yaml for ATAC or ChIP

### ATAC|RNA-seq design file

An ATAC-seq or RNA-seq design file should look something like this:

```bash
sample,r1,r2
rna,/path/to/fastq/rna_2.fastq.gz,/path/to/fastq/rna_1.fastq.gz
```

!!! Note
    The design file is a CSV file with the following columns:
      * `sample` - The sample name. Altering this will change the name of the output files so can be useful for renaming samples.
      * `r1` - The path to the read 1 fastq file
      * `r2` - The path to the read 2 fastq file


### ChIP-seq design file

A ChIP assay design file should look something like this:

```bash
sample_name,ip,control,ip_r1,ip_r2,control_r1,control_r2,norm_group
chip-rx,MLL,input,fastq/chip-rx_MLL_1.fastq.gz,fastq/chip-rx_MLL_2.fastq.gz,fastq/chip-rx_input_1.fastq.gz,fastq/chip-rx_input_2.fastq.gz,all
```

!!! Note
    The design file is a CSV file with the following columns:
      * `sample` - The sample name. Altering this will change the name of the output files so can be useful for renaming samples.
      * `ip_r1` - The path to the IP read 1 fastq file
      * `ip_r2` - The path to the IP read 2 fastq file
      * `control_r1` - The path to the control read 1 fastq file
      * `control_r2` - The path to the control read 2 fastq file
      * `ip` - The name of the IP sample
      * `control` - The name of the control sample


### RNA-seq design file

An RNA-seq design file should look something like this:

```bash
sample,r1,r2
rna,/path/to/fastq/rna_2.fastq.gz,/path/to/fastq/rna_1.fastq.gz
```

If you want to run DeSeq2, then you will need to add an additional column to the design file to indicate which samples are in the control group:

```bash
sample,r1,r2,deseq2
rna1,/path/to/fastq/rna1_2.fastq.gz,/path/to/fastq/rna1_1.fastq.gz,control
rna2,/path/to/fastq/rna2_2.fastq.gz,/path/to/fastq/rna2_1.fastq.gz,control
rna3,/path/to/fastq/rna3_2.fastq.gz,/path/to/fastq/rna3_1.fastq.gz,control
rna4,/path/to/fastq/rna4_2.fastq.gz,/path/to/fastq/rna4_1.fastq.gz,treated
rna5,/path/to/fastq/rna5_2.fastq.gz,/path/to/fastq/rna5_1.fastq.gz,treated
rna6,/path/to/fastq/rna6_2.fastq.gz,/path/to/fastq/rna6_1.fastq.gz,treated
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

screen -S NAME_OF_SESSION

# to exit screen session
  ctrl+a d 

# or 

tmux new -s NAME_OF_SESSION

# to detach from tmux session
  ctrl+b d

```

### Check you have activated the conda environment

```bash
conda activate seqnado
```

### Run the pipeline

The pipeline can be run using the following command:

```bash
seqnado [atac|chip|rna|snp] -c <number of cores> --preset [ss|ls] 
# additional options
--queue/-q [short|long] --scale-resource/-s <factor to multiply resources> 
```

An actual example would be:

```bash
seqnado rna -c 8 --preset ss -q short 
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
