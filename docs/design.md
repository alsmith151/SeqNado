[Back to Index](index.md)

After successful configuration of SeqNado run [Configuration](configuration.md)

## Design

The `seqnado design` command generates a design CSV file from FASTQ files for a specific assay. If no assay is provided, the tool operates in multiomics mode.

### Arguments

- **assay**: Specify the assay type. Options include:
  - `rna`
  - `atac`
  - `snp`
  - `chip`
  - `cat`
  - `meth`
  - `mcc`
  - `crispr`
  - If omitted, multiomics mode is used.
- **files**: Provide one or more FASTQ files.

### Options

- `--output`, `-o`: Specify the output CSV filename. Default: `metadata_{assay}.csv`.
- `--group-by`: Group samples by a regular expression or a column.
- `--auto-discover` / `--no-auto-discover`: Automatically search common folders if none are provided. Default: `auto-discover`.
- `--interactive` / `--no-interactive`: Interactively add missing columns using schema defaults. Default: `interactive`.
- `--accept-all-defaults`: Non-interactive mode that auto-adds only columns with schema defaults.
- `--verbose`, `-v`: Increase logging verbosity.
- `--help`: Display help information.

### Example Usage

#### Generate a Design CSV for RNA-seq
```bash
seqnado design rna fastqs/* 
```

#### Multiomics Mode
```bash
seqnado design 
```

## Next Steps

Once you have designed your experiment, proceed to set up the pipeline:

[Pipeline](pipeline.md)