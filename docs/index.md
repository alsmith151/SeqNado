# SeqNado

![Seqnado Logo](https://raw.githubusercontent.com/alsmith151/SeqNado/main/containers/seqnado.png)

SeqNado is a Snakemake-based workflow for ChIP-seq (with optional spike-in normalisation), ATAC-seq, RNA-seq, and short-read WGS SNP calling. It is designed to be modular, reproducible, and easy to deploy using Apptainer/Singularity containers. SeqNado provides end-to-end pipelines that take you from raw FASTQ files to publication-ready results with minimal setup.

## Quick Start

Follow the steps below to spin up a working pipeline instance. Consult the [installation](installation.md) page or other docs as needed.

### 1. Install SeqNado

```bash
mamba install -c bioconda seqnado
```

The installation guide covers alternative methods, optional dependencies, and troubleshooting tips.

### 2. Scaffold a working directory

```bash
seqnado-config [atac|chip|rna|snp]
```

This command creates a project directory and the configuration skeleton for the selected pipeline.

### 3. Stage input data

- Place or symlink raw FASTQ files into the generated `fastq/` folder.
- Copy or create any additional resources (e.g., spike-in references) expected by your config.

```bash
ln -s /path/to/fastq/*.fastq.gz /path/to/working-directory/fastq/
```

### 4. (Optional) Supply a design file

If your samples follow SeqNado naming conventions, the design table is generated automatically. Otherwise, run:

```bash
cd /path/to/working-directory
seqnado-design [atac|chip|rna|snp]
```

### 5. Launch the workflow

```bash
cd /path/to/working-directory
seqnado [atac|chip|rna|snp] -c <cores> --preset [ss|ls] \
  --queue/-q [short|long] --scale-resource/-s <factor>

# Example
seqnado rna -c 8 --preset ss
```

Use `--preset ss` to submit jobs to the cluster (recommended) or `--preset ls` for local execution. Adjust `--queue` and `--scale-resource` to fit your environment. Additional usage details, presets, and pipeline-specific parameters are documented in `docs/pipeline.md` and the CLI help (`seqnado --help`).
