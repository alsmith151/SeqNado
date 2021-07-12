# NGS Pipeline

Pipeline based on cgat-core (ruffus wrapper) to process ChIP-seq, ATAC-seq and RNA-seq data.

## Installation

1. Install conda if it has not been already using the [conda install instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent)

1b. (Optional) Install mamba. Mamba is a drop-in replacement for conda that is **much** faster than conda. All "conda" commands can be replaced with "mamba" for a   significant speedup.

```
conda install mamba
````

2. Clone the repository and enter it:

```
git clone https://github.com/alsmith151/ngs_pipeline.git
cd ngs_pipeline
```

3. Generate a new conda environment using the environment.yml file.

```
conda env create -f environment.yml
conda activate ngs
```

4. Install the pipeline:

```
pip install .
```

5. If you intend to use a cluster e.g. SLURM add the path to the DRMAA interface to your .bashrc:

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

2. Get and edit the pipeline configuration file.

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
# To use gedit
gedit config.yml

# To use nano
nano config.yml
```

3.  Copy or link fastq files into the working directory

The pipeline requires that fastq files are paired and in any of these formats:

For ChIP-seq:

Note that the underscore is needed to identify pairs for peak calling

* samplename1_Antibody_R1.fastq.gz
* samplename1_Antibody_R2.fastq.gz
* samplename1_Input_1.fastq
* samplename1_Input_2.fastq

For ATAC-seq: 

**Note the absence of an underscore separating the sample name.**

* sample-name-1_R1.fastq.gz
* sample-name-1_R2.fastq.gz
* sample-name-1_1.fastq
* sample-name-1_2.fastq

For RNA-seq:

* sample-name-1_R1.fastq.gz
* sample-name-1_R2.fastq.gz
* sample-name-1_1.fastq
* sample-name-1_2.fastq

All FASTQ files present in the directory will be processed by the pipeline in parallel and
original FASTQ files will not be modified. If new FASTQ files are added to a pre-run pipeline,
only the new files will be processed.

Copy:
```
cp PATH_TO_FASTQ/example_R1.fastq.gz
```

Symlink:
```
# Be sure to use the absolute path for symlinks
ln -s /ABSOLUTE_PATH_TO_FASTQ/example_R1.fastq.gz
```

4. Running the pipeline

After copying/linking FASTQ files into the working directory and configuring the copy of
config.yml in the working directory for the current experiment, the pipeline can be run with:

```
ngs-pipeline atac # ATAC-seq samples
ngs-pipeline chip # ChIP-seq/ChIPMentation
ngs-pipeline rna # RNA-seq - Not fully tested 
```

There are several options to visualise which tasks will be performed by the pipeline
before running. 

The tasks to be performed can be examined with:
```    
# Shows the tasks to be performed
ngs-pipeline atac show 

# Plots a directed graph using graphviz
ngs-pipeline atac plot
```

If you are happy with the tasks to be performed, the full pipeline run can be launched with:

```
# If using all default settings and using a cluster
ngs-pipeline atac make

# Higher verbosity
ngs-pipeline atac make -v 5

# If not using a cluster, run in local mode.
ngs-pipeline atac make --local -p 4

# Avoiding network disconnections
nohup ngs-pipeline atac make &
```

See [cgat-core Read the Docs](https://cgat-core.readthedocs.io/en/latest/getting_started/Examples.html) for additional
information.







