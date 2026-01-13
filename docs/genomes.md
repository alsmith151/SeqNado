[← Back to main page](index.md)

# Genome Setup

After successful initialization of SeqNado ([Initialisation](initialisation.md)).

The `seqnado genomes` command allows you to manage genome configurations, including listing, editing, building, and generating configurations for `fastq-screen`.

## Examples

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

For all subcommands, arguments, and flags, see the CLI reference: [seqnado genomes](cli.md#cli-seqnado-genomes).

## Next Steps

Once you have configured your genomes, proceed to Configure a SeqNado run:

[Configuration](configuration.md)

