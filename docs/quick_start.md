[Back to Index](index.md)

# Quick Start

Get started with SeqNado in just a few steps.

SeqNado can be run for any of the following assay types, as well as in multiomics mode:

- ATAC-seq: `atac`
- ChIP-seq: `chip`
- CUT&Tag: `cat`
- RNA-seq: `rna`
- Methylation: `meth`
- SNP analysis: `snp`
- MCC: `mcc`
- CRISPR analysis: `crispr`

## Example Workflow
### 1. Install SeqNado

The fastest method to install SeqNado is from bioconda via mamba.

```bash
mamba create -n seqnado -c bioconda seqnado
mamba activate seqnado
```

### 2. Initialise SeqNado

The `seqnado init` command initializes the SeqNado user environment. This step ensures that the necessary configuration files and dependencies are set up for the package to function correctly.

#### Usage
```bash
seqnado init [OPTIONS]
```

#### Options
- `--preset, --no-preset`: Use packaged preset genomes instead of the editable template (default: disabled).
- `--dry-run, --no-dry-run`: Show actions without writing files or running scripts (default: disabled).
- `--verbose, -v`: Increase logging verbosity.

#### Actions Performed
- Logs the current Conda environment if active (optional).
- Runs the packaged Apptainer/Singularity initialization if `apptainer` is available on the system PATH.
- Ensures the `~/.config/seqnado/genome_config.json` file exists, either as a template or using a preset.

#### Example
Initialize SeqNado with default settings:
```bash
seqnado init
```

### 3. Set up genome references for SeqNado

The `seqnado genomes` command manages genome configurations, including listing, editing, building, or generating `fastq-screen` configurations.

#### Usage
```bash
seqnado genomes [OPTIONS] SUBCOMMAND ASSAY
```

#### Subcommands
- **list**: List available genome configurations.
- **edit**: Edit an existing genome configuration.
- **build**: Build a new genome configuration from a FASTA file.
- **fastqscreen**: Generate a `fastq-screen` configuration file.

#### Arguments
- **SUBCOMMAND**: The operation to perform (e.g., `list`, `edit`, `build`, `fastqscreen`).
- **ASSAY**: Assay type. Options include `rna`, `atac`, `snp`, `chip`, `cat`, `meth`, `mcc`, `crispr` (default: `atac`).

#### Options
- `--fasta, -f`: Input FASTA file (required for `build`).
- `--name, -n`: Genome name (prefix) for the built genome.
- `--outdir, -o`: Output directory for genome builds (default: `/ceph/project/milne_group/cchahrou/software/SeqNado/genome_build`).
- `--screen, -s`: Output path for `fastq-screen` configuration files (used in the `fastqscreen` subcommand).
- `--threads, -t`: Number of threads for Bowtie2 (default: 8).
- `--no-contaminants`: Exclude contaminant databases in `fastq-screen` configurations.
- `--contaminant-path`: Path to contaminant reference files.
- `--verbose, -v`: Increase logging verbosity.

#### Example
Build a genome configuration for RNA-seq:
```bash
seqnado genomes build rna --fasta hg38.fasta --name hg38 --outdir /path/to/output
```

### 4. Configure a SeqNado run

The `seqnado config` command builds a workflow configuration YAML for the selected assay. If no assay is provided, the command operates in multiomics mode.

#### Usage
```bash
seqnado config [OPTIONS] [ASSAY]
```

#### Arguments
- **ASSAY**: Assay type. Options include `rna`, `atac`, `snp`, `chip`, `cat`, `meth`, `mcc`, `crispr`. If omitted, multiomics mode is used.



#### Options
- `--make-dirs, --no-make-dirs`: Create or skip creating the output project directory or FASTQ subdirectory (default: create).
- `--render-options, --no-render-options`: Render all options, even if not used by the workflow (default: disabled).
- `--output, -o`: Specify the explicit path for the rendered configuration file.
- `--interactive, --no-interactive`: Enable or disable interactive prompts for configuration values (default: enabled).
- `--verbose, -v`: Increase logging verbosity.

#### Example
Generate a configuration file:
```bash
seqnado config 
```

You can edit the generated YAML file to customize the workflow for your specific needs.

### 5. SeqNado design

The `seqnado design` command generates a metadata design CSV from FASTQ files for a specific assay. If no assay is provided, the command operates in multiomics mode. The generated CSV outlines the structure of the experiment, including sample names, conditions, and other relevant metadata.

#### Usage
```bash
seqnado design [OPTIONS] [ASSAY] [FASTQ ...]
```

#### Arguments
- **ASSAY**: Assay type. Options include `rna`, `atac`, `snp`, `chip`, `cat`, `meth`, `mcc`, `crispr`. If omitted, multiomics mode is used.
- **FASTQ**: One or more FASTQ files to include in the design.

#### Options
- `--output, -o`: Specify the output CSV filename (default: `metadata_{assay}.csv`).
- `--group-by`: Group samples by a regular expression or a column.
- `--auto-discover`: Automatically search common folders for FASTQ files if none are provided (default: enabled).
- `--interactive`: Interactively add missing columns using schema defaults (default: enabled).
- `--accept-all-defaults`: Non-interactive mode; auto-add only columns with schema defaults.
- `--verbose, -v`: Increase logging verbosity.

#### Example
Generate a design CSV for ATAC-seq:
```bash
seqnado design
```

The generated CSV can be reviewed and edited to ensure all experimental details are correctly specified.

### 6. Run SeqNado pipeline

The `seqnado pipeline` command runs the data processing pipeline for the specified assay. It uses Snakemake under the hood to manage the workflow.

#### Usage
```bash
seqnado pipeline [OPTIONS] [ASSAY]
```

#### Arguments
- **ASSAY**: Assay type. Required for single-assay workflows, optional for multiomics mode.

#### Options
- `--configfile`: Path to a SeqNado configuration YAML file (default: `config_<ASSAY>.yaml`).
- `--preset`: Snakemake job profile preset. Options include:
  - `lc`: Local cluster
  - `le`: Local execution (default)
  - `ls`: Local single-threaded
  - `ss`: Slurm scheduler
  - `t`: Test mode
- `--clean-symlinks, --no-clean-symlinks`: Remove symlinks created by previous runs (default: disabled).
- `--scale-resources, -s`: Scale memory and time resources (default: 1.0).
- `--verbose, -v`: Increase logging verbosity.
- `--queue, -q`: Specify the Slurm queue/partition for the `ss` preset (default: `short`).
- `--print-cmd`: Print the Snakemake command before running it.

#### Example
Run the pipeline for ATAC-seq with 8 cores and the `ls` preset:
```bash
seqnado pipeline atac --preset ls --scale-resources 1.5
```
