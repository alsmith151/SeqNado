# NGS Pipeline

Pipeline based on cgat-core (ruffus wrapper) to process ChIP-seq and ATAC-seq data.

## Installation

1. Install conda if it has not been already using the [conda install instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent)

1b. (Optional) Install mamba. Mamba is a drop-in replacement for conda that is **much** faster than conda. All "conda" commands can be replaced with "mamba" for a significant speedup.

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
export DRMAA_LIBRARY_PATH=/usr/lib64/libdrmaa.so
```






