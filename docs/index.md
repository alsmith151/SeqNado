![Seqnado Logo](https://raw.githubusercontent.com/alsmith151/SeqNado/master/seqnado_logo.jpeg)

Pipeline based on snakemake to process ChIP-seq (with optional spike-in based normalisation), ATAC-seq, RNA-seq and short read WGS data for SNP calling. The defaults are optimised for the Milne group directory on the CCB cluster, but can be easily modified for other groups and clusters.


## Quick Start

### Installation

The [installation](installation.md) page has detailed instructions for installing SeqNado. For a very quick start, run the following command:

```bash
bash <(curl -s https://raw.githubusercontent.com/alsmith151/SeqNado/master/install_seqnado.sh)
```

### Set up the pipeline

#### Config file

Generate the config file and a working directory for the pipeline:

```bash
seqnado-config [atac|chip|rna|snp]
```

#### Design (optional)

The design will be generated automatically if not provided and the samples follow the correct naming convention.

```bash
seqnado-design [atac|chip|rna|snp]
```

#### Ensure files are in the correct location

Ensure that the fastq files, and design are in the correct location:

```bash
# Fastq files
ln -s /path/to/fastq_files/ /path/to/working-directory/made-by-seqnado-config/

# Design
mv /path/to/design.csv /path/to/working-directory/made-by-seqnado-config/
```

### Run the pipeline

```bash
cd /path/to/working-directory/made-by-seqnado-config/
seqnado [atac|chip|rna|snp] -c <number of cores> --preset [ss|ls] # ss = use cluster, ls = use local (not recommended)
```
