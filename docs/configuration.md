[Back to Index](index.md)

# Configuration

After successful genome configuration of SeqNado [Genomes](genomes.md)


## Using `seqnado config`

The `seqnado config` command builds a workflow configuration YAML file for the selected assay. If no assay is provided, the tool operates in multiomics mode.

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

### Options

- `--make-dirs` / `--no-make-dirs`: Create or skip creating the output project directory or fastq subdirectory. Default: `make-dirs`.
- `--render-options` / `--no-render-options`: Render all options, even if not used by the workflow. Default: `no-render-options`.
- `--output`, `-o`: Specify the explicit path for the rendered configuration file.
- `--verbose`, `-v`: Increase logging verbosity.
- `--interactive` / `--no-interactive`: Interactively prompt for configuration values. Default: `interactive`.
- `--help`: Display help information.

### Example Usage

#### Build a Configuration for RNA-seq
```bash
seqnado config rna --output rna_config.yaml
```

#### Multiomics Mode
```bash
seqnado config --make-dirs --interactive
```

### FASTQ Files

After generating the configuration and project directory using `seqnado config`, you need to link your FASTQ files into the `fastqs` directory. This ensures that the pipeline can locate and process your input data.

#### Symlinking FASTQ Files

Use the following command to create symbolic links for your FASTQ files:

```bash
ln -s /path/to/your/fastq/files/* <project_directory>/fastqs/
```

Replace `/path/to/your/fastq/files/` with the directory containing your FASTQ files and `<project_directory>` with the path to the project directory created by `seqnado config`.

#### Example

If your FASTQ files are located in `/data/fastq/` and your project directory is `rna_project`, run:

```bash
ln -s /data/fastq/* rna_project/fastqs/
```

This will create symbolic links to all FASTQ files in the `fastqs` directory of your project.

#### Safe Naming Strategies for FASTQ Files

To ensure compatibility with the pipeline and avoid errors, use consistent and descriptive naming conventions for your FASTQ files. Below are some examples of safe naming strategies:


- **ATAC-seq**:
  ```
  sample-name-rep1_R1.fastq.gz
  sample-name-rep1_R2.fastq.gz
  ```

- **ChIP-seq**:
  ```
  sample-name-rep1_Antibody_R1.fastq.gz
  sample-name-rep1_Antibody_R2.fastq.gz
  sample-name-rep2_Input_R1.fastq.gz
  sample-name-rep2_Input_R2.fastq.gz
  ```
  - `Antibody`: Name of the antibody used for ChIP.
  - `Input`: Control sample.

- **RNA-seq**:
  ```
  sample-name-rep1_R1.fastq.gz
  sample-name-rep1_R2.fastq.gz
  sample-name-rep2_R1.fastq.gz
  sample-name-rep2_R2.fastq.gz
  ```
  - `sample-name`: Unique identifier for the sample.
  - `rep1`, `rep2`: Biological or technical replicate number.
  - `R1`, `R2`: Read pair (forward and reverse).

Using these naming conventions ensures that the pipeline can correctly parse and process your data.

### Third Party Tools

SeqNado integrates with several third-party tools to enhance its functionality. 
The configuration yaml file can be edited to set the parameters for use.
Below is a categorized list of supported tools:

#### Alignment Tools
- **bowtie2**: For aligning sequencing reads to a reference genome.
- **star**: A fast aligner for RNA-seq data.

#### Processing Tools
- **samtools**: Utilities for manipulating alignments in the SAM/BAM format.
- **picard**: A set of Java command-line tools for manipulating high-throughput sequencing data.

#### Trimming Tools
- **cutadapt**: Removes adapter sequences from high-throughput sequencing reads.
- **trimgalore**: A wrapper tool around Cutadapt and FastQC for quality control and adapter trimming.

#### Peak Calling Tools
- **macs**: Model-based analysis of ChIP-seq data.
- **lanceotron**: A tool for identifying peaks in ChIP-seq data.
- **lanceotronmcc**: Specialized for MCC data.
- **seacr**: A peak caller for sparse data.

#### Analysis Tools
- **deeptools**: For the analysis and visualization of high-throughput sequencing data.
- **homer**: A suite of tools for motif discovery and next-gen sequencing analysis.
- **bamnado**: A tool for BAM file manipulation.
- **subread**: A tool for read alignment and feature quantification.

#### Specialized Tools
- **methyldackel**: For processing bisulfite sequencing data.
- **bcftools**: Utilities for variant calling and manipulating VCF/BCF files.

For more details on configuring these tools, refer to the [Third Party Tools Configuration Guide](third_party_tools_config.md).

## Next Steps

Once you have configured SeqNado, proceed to design your workflows:

[Design](design.md)