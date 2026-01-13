<!-- Header and badges -->
# SeqNado

[![Documentation](https://github.com/alsmith151/SeqNado/actions/workflows/build_docs.yml/badge.svg)](https://github.com/alsmith151/SeqNado/actions/workflows/build_docs.yml)
[![Bioconda](https://anaconda.org/bioconda/seqnado/badges/version.svg)](https://anaconda.org/bioconda/seqnado)
[![Bioconda Updated](https://anaconda.org/bioconda/seqnado/badges/latest_release_date.svg)](https://anaconda.org/bioconda/seqnado)
[![PyPI Downloads](https://static.pepy.tech/badge/seqnado)](https://pepy.tech/projects/seqnado)

<p align="center">
  <img src="https://raw.githubusercontent.com/alsmith151/SeqNado/main/containers/pipeline/seqnado.png" alt="SeqNado logo" />
</p>

*A powerful, unified bioinformatics toolkit for ATAC-seq, ChIP-seq, CUT&Tag, RNA-seq, SNP analysis, Methylation, CRISPR screens, and Micro-Capture-C.*

## Table of Contents

- [Key Features](#key-features)
- [Supported Assays](#supported-assays)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Documentation](#documentation)
- [License](#license)

---

Modular, reproducible, and container-ready pipelines powered by Snakemake that take you from raw data to publication-ready results.

## Key Features

- **Comprehensive Assay Support**: Single framework for multiple sequencing assays
- **User-Friendly CLI**: Intuitive command-line interface that guides you through setup and execution
- **Multiomics Support**: Analyze and integrate data from multiple sequencing assays in a single workflow
- **Snakemake-Powered**: Modular workflows with automatic parallelization and resource management
- **Container-Ready**: Fully containerized pipelines using Apptainer/Singularity for reproducibility
- **HPC-Optimized**: Seamless integration with SLURM and local execution modes
- **Advanced Analysis**: 
  - Spike-in normalization for ChIP-seq, ATAC-seq, and RNA-seq
  - Automated differential expression with DESeq2
  - UCSC genome browser hub generation
  - Peak calling with MACS2, SEACR, and LanceOtron
  - Comprehensive QC with MultiQC reports
- **Flexible Configuration**: Interactive CLI for setup, or scriptable non-interactive mode
- **Machine Learning Ready**: Tools for preparing datasets for ML applications

## Supported Assays

- **ATAC-seq** (`atac`) - Chromatin accessibility profiling with TSS enrichment and fragment analysis
- **ChIP-seq** (`chip`) - Protein-DNA interaction mapping with spike-in support
- **CUT&Tag** (`cat`) - Low-input epigenomic profiling optimized for sparse signals
- **RNA-seq** (`rna`) - Transcriptome analysis with automated DESeq2 differential expression
- **SNP Analysis** (`snp`) - Variant detection and genotyping workflows
- **Methylation** (`meth`) - Bisulfite/TAPS sequencing for DNA methylation analysis
- **CRISPR Screens** (`crispr`) - Guide-level quantification and screen statistics
- **Micro-Capture-C** (`mcc`) - Micro-Capture-C chromatin conformation capture analysis
- **Multiomics** - Integrated analysis in parallel across multiple assay types

## Installation

### Via Mamba (Recommended)

Install from the Bioconda channel:

```bash
mamba create -n seqnado -c bioconda seqnado
mamba activate seqnado
```

### Via Pip

Alternatively, install from PyPI:

```bash
pip install seqnado
```

### Initialize SeqNado

After installation, initialize your SeqNado environment:

```bash
seqnado init
```

This sets up genome configurations, Apptainer containers, and Snakemake profiles.

## Quick Start

### 1. Set Up Genome References

List available genomes or build a new one:

```bash
# List available genomes
seqnado genomes list atac

# Build a custom genome
seqnado genomes build rna --fasta hg38.fasta --name hg38 --outdir /path/to/genomes
```

### 2. Create Project Configuration

Generate a configuration file for your assay:

```bash
seqnado config atac
```

This creates a project directory with the structure:
```
YYYY-MM-DD_ASSAY_project/
├── config_atac.yaml
└── fastqs/
```

### 3. Add FASTQ Files

Link your FASTQ files into the project:

```bash
ln -s /path/to/fastq/*.fastq.gz YYYY-MM-DD_ATAC_project/fastqs/
```

### 4. Generate Sample Metadata

Create the experimental design:

```bash
seqnado design atac
```

This generates a `metadata_atac.csv` file that you can edit to specify sample groupings, conditions, and controls.

### 5. Run the Pipeline

Execute the workflow:

```bash
# Local execution
seqnado pipeline atac --preset le

# HPC cluster with SLURM
seqnado pipeline atac --preset ss --queue short

# Multiomics mode (if multiple config files present)
seqnado pipeline
```

### Common Pipeline Options

- `--preset le` - Local execution (default)
- `--preset ss` - SLURM scheduler for HPC
- `--preset lc` - Local execution with conda environment
- `--scale-resources 1.5` - Scale memory/time by 1.5x
- `--queue short` - Specify SLURM partition
- `-n` - Dry run (passed to Snakemake)
- `--unlock` - Unlock working directory (passed to Snakemake)

## Documentation

For comprehensive guides and API documentation, visit:

**📚 [SeqNado Documentation](https://alsmith151.github.io/SeqNado/)**

### Key Topics

- [Installation Guide](https://alsmith151.github.io/SeqNado/installation/)
- [Genome Setup](https://alsmith151.github.io/SeqNado/genomes/)
- [Configuration](https://alsmith151.github.io/SeqNado/configuration/)
- [Design Files](https://alsmith151.github.io/SeqNado/design/)
- [Pipeline Details](https://alsmith151.github.io/SeqNado/pipeline/)
- [Outputs](https://alsmith151.github.io/SeqNado/outputs/)
- [CLI Reference](https://alsmith151.github.io/SeqNado/cli/)
- [HPC Cluster Configuration](https://alsmith151.github.io/SeqNado/cluster_config/)

## License

This project is licensed under [GPL3](LICENSE).
