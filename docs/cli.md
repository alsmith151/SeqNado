[← Back to main page](index.md)

# CLI Reference

This page provides a complete reference of all available SeqNado commands and options.

For quick examples and typical usage patterns, see the [Quick Start Guide](quick_start.md).

---

# `seqnado`

**SeqNado CLI**

Initialize your environment, build configs, create design files, and run pipelines.
Use --help on any subcommand for details.

**Usage**:

```console
$ seqnado [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `-v, --version`: Show version and exit.
* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `init`: Initialize SeqNado user environment.
* `genomes`: Manage genome configurations (list, edit, build, fastqscreen)
* `config`: Build a workflow configuration YAML
* `design`: Generate a SeqNado design CSV from FASTQ files
* `pipeline`: Run the data processing pipeline

## `seqnado init` {#cli-seqnado-init}

Initialize SeqNado user environment.

- Logs the current Conda environment if active (optional).
- Runs packaged Apptainer/Singularity init (if `apptainer` on PATH).
- Ensures ~/.config/seqnado/genome_config.json exists (template or preset).

**Usage**:

```console
$ seqnado init [OPTIONS]
```

**Options**:

* `--preset / --no-preset`: Use packaged preset genomes instead of the editable template.  [default: no-preset]
* `--dry-run / --no-dry-run`: Show actions without writing files or running scripts.  [default: no-dry-run]
* `-v, --verbose`: Increase logging verbosity.
* `--help`: Show this message and exit.

## `seqnado genomes` {#cli-seqnado-genomes}

Manage genome configurations (list, edit, build, or generate fastq-screen config)

**Usage**:

```console
$ seqnado genomes [OPTIONS] SUBCOMMAND [ASSAY]
```

**Arguments**:

- **SUBCOMMAND**: `list` | `edit` | `build` | `fastqscreen` (required)
- **ASSAY**: `rna` | `atac` | `snp` | `chip` | `cat` | `meth` | `mcc` | `crispr` | `multiomics` (default: `atac`)

**Options**:

* `-f, --fasta PATH`: Input FASTA (required for build)
* `-n, --name TEXT`: Genome name (prefix) for built genome
* `-o, --outdir PATH`: Output directory for build
* `-s, --screen PATH`: Output path for fastqscreen config (fastqscreen subcommand)
* `-t, --threads INTEGER`: Number of threads for Bowtie2 (fastqscreen subcommand)  [default: 8]
* `--no-contaminants`: Exclude contaminant databases (fastqscreen subcommand)
* `--contaminant-path PATH`: Path to contaminant reference files (fastqscreen subcommand)
* `-v, --verbose`: Increase logging verbosity
* `--help`: Show this message and exit.

## `seqnado config` {#cli-seqnado-config}

Build a workflow configuration YAML for the selected ASSAY. If no assay is provided, multiomics mode is used.

**Usage**:

```console
$ seqnado config [OPTIONS] [ASSAY]
```

**Arguments**:

- **ASSAY**: `rna` | `atac` | `snp` | `chip` | `cat` | `meth` | `mcc` | `crispr` | `multiomics`; omitted → multiomics mode

**Options**:

* `--make-dirs / --no-make-dirs`: Create/don't create the output project directory or fastq subdir.  [default: make-dirs]
* `--render-options / --no-render-options`: Render all options (even if not used by the workflow).  [default: no-render-options]
* `-o, --output PATH`: Explicit path for the rendered config file.
* `-v, --verbose`: Increase logging verbosity.
* `--interactive / --no-interactive`: Interactively prompt for config values. Non-interactive mode only works for single assay configs (except MCC and multiomics).  [default: interactive]
* `--help`: Show this message and exit.

## `seqnado design` {#cli-seqnado-design}

Generate a SeqNado design CSV from FASTQ files for ASSAY. If no assay is provided, multiomics mode is used.

**Usage**:

```console
$ seqnado design [OPTIONS] [ASSAY] [FASTQ ...]
```

**Arguments**:

- **ASSAY**: `rna` | `atac` | `snp` | `chip` | `cat` | `meth` | `mcc` | `crispr` | `multiomics`; omitted → multiomics mode
- **FASTQ**: one or more FASTQ files

**Options**:

* `-o, --output PATH`: Output CSV filename (default: metadata_{assay}.csv).
* `--ip-to-control TEXT`: List of antibody,control pairings for IP assays (e.g. ChIP). Format: 'antibody1:control1,antibody2:control2'. If provided will assign a control with a specified name to that ip in the metadata. If not provided, controls will be assigned based on a best-effort matching of sample names.
* `--group-by`: Group samples by a regular expression or a column.
* `--auto-discover / --no-auto-discover`: Search common folders if none provided.  [default: auto-discover]
* `--interactive / --no-interactive`: Interactively offer to add missing columns using schema defaults.  [default: interactive]
* `--accept-all-defaults`: Non-interactive: auto-add only columns that have a schema default.
* `--deseq2-pattern TEXT`: Regex pattern to extract DESeq2 groups from sample names. First capture group will be used. Example: r'-(\w+)-rep' for 'sample-GROUP-rep1'
* `-v, --verbose`: Increase logging verbosity.
* `--help`: Show this message and exit.

## `seqnado pipeline` {#cli-seqnado-pipeline}

Run the data processing pipeline for ASSAY (Snakemake under the hood). Any additional arguments are passed to Snakemake (e.g., `seqnado pipeline rna -n` for dry-run, `--unlock`, etc.).

**Usage**:

```console
$ seqnado pipeline [OPTIONS] [ASSAY]
```

**Arguments**:

- **ASSAY**: required for single-assay; optional in multiomics mode

**Options**:

* `--configfile PATH`: Path to a SeqNado config YAML (default: config_<ASSAY>.yaml).
* `--version`: Print SeqNado version and exit.
* `--preset [lc|le|ls|ss|t]`: Snakemake job profile preset.  [default: le]
* `--clean-symlinks / --no-clean-symlinks`: Remove symlinks created by previous runs.  [default: no-clean-symlinks]
* `-s, --scale-resources FLOAT`: Scale memory/time (env: SCALE_RESOURCES).  [default: 1.0]
* `-v, --verbose`: Increase logging verbosity.
* `-q, --queue TEXT`: Slurm queue/partition for the `ss` preset.  [default: short]
* `--print-cmd`: Print the Snakemake command before running it.
* `--help`: Show this message and exit.
