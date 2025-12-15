[Back to Index](index.md)

# Managing Genome Configurations

After successful initialization of SeqNado [Initialisation](initialisation.md)

The `seqnado genomes` command allows you to manage genome configurations, including listing, editing, building, and generating configurations for `fastq-screen`.

## Subcommands

- **list**: Display available genome configurations.
- **edit**: Modify an existing genome configuration.
- **build**: Create a new genome configuration using a FASTA file.
- **fastqscreen**: Generate a configuration file for `fastq-screen`.

## Arguments

- **subcommand**: Required. Choose one of the subcommands: `list`, `edit`, `build`, or `fastqscreen`.
- **assay**: Optional. Specify the assay type. Options include:
  - `rna`
  - `atac` (default)
  - `snp`
  - `chip`
  - `cat`
  - `meth`
  - `mcc`
  - `crispr`

## Options

- `--fasta`, `-f`: Path to the input FASTA file (required for the `build` subcommand).
- `--name`, `-n`: Genome name (prefix) for the built genome.
- `--outdir`, `-o`: Output directory for the genome build. Default: `/ceph/project/milne_group/cchahrou/software/SeqNado/genome_build`.
- `--screen`, `-s`: Output path for the `fastq-screen` configuration file (used with the `fastqscreen` subcommand).
- `--threads`, `-t`: Number of threads for Bowtie2 (used with the `fastqscreen` subcommand). Default: 8.
- `--no-contaminants`: Exclude contaminant databases (used with the `fastqscreen` subcommand).
- `--contaminant-path`: Path to contaminant reference files (used with the `fastqscreen` subcommand).
- `--verbose`, `-v`: Increase logging verbosity.
- `--help`: Display help information.

## Example Usage

### List Available Genomes
```bash
seqnado genomes list
```

### Build a New Genome
```bash
seqnado genomes build --fasta /path/to/genome.fasta --name my_genome --outdir /path/to/output
```

### Generate a `fastq-screen` Configuration
```bash
seqnado genomes fastqscreen --screen /path/to/fastqscreen.conf --threads 16 --no-contaminants
```

## Next Steps

Once you have configured your genomes, proceed to Configure a SeqNado run:

[Configuration](configuration.md)

