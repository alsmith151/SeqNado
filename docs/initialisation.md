[Back to Index](index.md)

# Initializing SeqNado

After successful installation of SeqNado [Installation](installation.md)

## Initialize the Workflow
Run the following command to initialize the SeqNado workflow in your project directory:

```bash
seqnado init
```

This will create the necessary configuration files and directory structure for your environment.

## Additional Initialization Options

The `seqnado init` command provides several options to customize the initialization process:

- `--preset` / `--no-preset`: Use packaged preset genomes instead of the editable template. Default is `--no-preset`.
- `--dry-run` / `--no-dry-run`: Show actions without writing files or running scripts. Default is `--no-dry-run`.
- `--verbose` / `-v`: Increase logging verbosity for detailed output.
- `--help`: Display help information for the `seqnado init` command.

### Example Usage

To initialize SeqNado with preset genomes and verbose logging:

```bash
seqnado init --preset --verbose
```

To perform a dry run and preview the initialization steps:

```bash
seqnado init --dry-run
```

## Next Steps

Once SeqNado has been Initialized, continue with:

[Genomes](genomes.md)