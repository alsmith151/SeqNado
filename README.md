<!-- Header and badges -->
# SeqNado Pipeline

[![Documentation](https://github.com/alsmith151/SeqNado/actions/workflows/build_docs.yml/badge.svg)](https://github.com/alsmith151/SeqNado/actions/workflows/build_docs.yml)
[![Bioconda](https://anaconda.org/bioconda/seqnado/badges/version.svg)](https://anaconda.org/bioconda/seqnado)
[![Bioconda Updated](https://anaconda.org/bioconda/seqnado/badges/latest_release_date.svg)](https://anaconda.org/bioconda/seqnado)
[![PyPI Downloads](https://static.pepy.tech/badge/seqnado)](https://pepy.tech/projects/seqnado)

![SeqNado logo](https://raw.githubusercontent.com/alsmith151/SeqNado/master/seqnado_logo.jpeg)

*A unified, user-friendly collection of workflows for ATAC-seq, ChIP-seq, CUT&RUN/TAG, RNA-seq, WGS, Methylation (Bisulfite/TAPS), CRISPR screens & Micro-Capture-C.*

## Table of Contents

- [Key Features](#key-features)
- [Installation](#installation)
- [Quickstart](#quickstart)
- [License](#license)

---

*A unified, user-friendly collection of Snakemake-based workflows for ATAC-seq, ChIP-seq, CUT&RUN/TAG, RNA-seq, whole-genome sequencing (WGS), methylation (Bisulfite/TAPS), CRISPR screens, and Micro-Capture-C.*

---

Empower your genomics research with modular, reproducible, and container-ready pipelines that take you from raw data to publication-ready results.

See the [SeqNado documentation](https://alsmith151.github.io/SeqNado/) for more information.

## Key Features

- Modular, Snakemake-based workflows for multiple assays
- Container-ready pipelines using Apptainer/Singularity
- Configurable and reproducible analysis environments
- Support for ATAC-seq, ChIP-seq, RNA-seq, WGS, methylation, CRISPR screens, and M-Capture-C

## Installation

Install via mamba (Bioconda channel):  

```bash
mamba install -c bioconda seqnado
```


## Quickstart

Generate a config and working directory for your assay:  

```bash
seqnado-config [atac|chip|rna|snp|meth|crispr|mcc]
```

Link your FASTQ files into the new directory:  

```bash
ln -s /path/to/fastq/* <workdir>/fastq/
```

Run the pipeline (example uses RNA-seq with 8 cores on the cluster preset):  

```bash
cd <workdir>
seqnado rna -c 8 --preset ss
```

For full details and advanced options, see the [SeqNado documentation](https://alsmith151.github.io/SeqNado/docs/index.html#quick-start).

## License

This project is licensed under the [MIT License](LICENSE).
