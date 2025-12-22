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
- `--ip-to-control`: List of antibody,control pairings for IP assays (e.g. ChIP). Format: 'antibody1:control1,antibody2:control2'.
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

#### Generate a Design CSV for ChIP-seq with explicit control pairing

For ChIP-seq experiments, the design CSV requires both IP FASTQ files and optionally control FASTQ files. The tool can infer these relationships based on file naming conventions.

#### Simple Case

For simple cases with a single control or when no control is needed, for example:

* SAMPLE1_H3K27ac_R1.fastq.gz
* SAMPLE1_H3K27ac_R2.fastq.gz
* SAMPLE1_Menin_R1.fastq.gz
* SAMPLE1_Menin_R2.fastq.gz
* SAMPLE1_input_R1.fastq.gz
* SAMPLE1_input_R2.fastq.gz
* SAMPLE_2_H3K27ac_R1.fastq.gz
* SAMPLE_2_H3K27ac_R2.fastq.gz

The command would be:

```bash
seqnado design chip fastqs/* 
```

The control will either be left blank if no appropriate files are in the directory or a single control sharing the same sample ID will be broadcast to all IP samples sharing that sample ID. e.g.:

| assay | sample_id | ip      | control | r1                           | r2                           | r1_control                | r2_control                | scaling_group |
|-------|-----------|---------|---------|------------------------------|------------------------------|---------------------------|---------------------------|---------------|
| ChIP  | SAMPLE1   | H3K27ac | input   | SAMPLE1_H3K27ac_R1.fastq.gz  | SAMPLE1_H3K27ac_R2.fastq.gz  | SAMPLE1_input_R1.fastq.gz | SAMPLE1_input_R2.fastq.gz | default       |
| ChIP  | SAMPLE1   | Menin   | input   | SAMPLE1_Menin_R1.fastq.gz    | SAMPLE1_Menin_R2.fastq.gz    | SAMPLE1_input_R1.fastq.gz | SAMPLE1_input_R2.fastq.gz | default       |
| ChIP  | SAMPLE_2  | H3K27ac |         | SAMPLE_2_H3K27ac_R1.fastq.gz | SAMPLE_2_H3K27ac_R2.fastq.gz |                           |                           | default       |



#### Complex Case with Multiple Controls and ambiguity in pairing

If there are multiple controls, specify which control corresponds to each IP using the `--ip-to-control` option. For example:

 We want the single fixed control `sf-input` to be used for the `H3K27ac` IP, and the double fixed `df-input` control to be used for the `Menin` IP. The FASTQ files are as follows:

  * SAMPLE1_H3K27ac_R1.fastq.gz
  * SAMPLE1_H3K27ac_R2.fastq.gz 
  * SAMPLE1_sf-input_R1.fastq.gz
  * SAMPLE1_sf-input_R2.fastq.gz
  * SAMPLE1_Menin_R1.fastq.gz
  * SAMPLE1_Menin_R2.fastq.gz
  * SAMPLE1_df-input_R1.fastq.gz  
  * SAMPLE1_df-input_R2.fastq.gz

The command would be:

```bash
seqnado design chip fastqs/* --ip-to-control "H3K4me3:sf-input,Menin:df-input"
```

This will generate a design CSV file with the appropriate control pairings. e.g.,

| assay | sample_id | ip      | control  | r1                          | r2                          | r1_control                   | r2_control                   | scaling_group |
|-------|-----------|---------|----------|-----------------------------|-----------------------------|------------------------------|------------------------------|---------------|
| ChIP  | SAMPLE1   | H3K27ac | sf-input | SAMPLE1_H3K27ac_R1.fastq.gz | SAMPLE1_H3K27ac_R2.fastq.gz | SAMPLE1_sf-input_R1.fastq.gz | SAMPLE1_sf-input_R2.fastq.gz | default       |
| ChIP  | SAMPLE1   | Menin   | df-input | SAMPLE1_Menin_R1.fastq.gz   | SAMPLE1_Menin_R2.fastq.gz   | SAMPLE1_df-input_R1.fastq.gz | SAMPLE1_df-input_R2.fastq.gz | default       |



#### Multiomics Mode
```bash
seqnado design 
```

## Next Steps

Once you have designed your experiment, proceed to set up the pipeline:

[Pipeline](pipeline.md)