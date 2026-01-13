[← Back to main page](index.md)

# Running SeqNado on HPC Clusters

SeqNado is designed to run efficiently on high-performance computing (HPC) clusters using job schedulers. This guide covers setting up and running SeqNado on SLURM-based clusters, which is the recommended and most thoroughly tested configuration.

## Quick Start: SLURM Clusters

To run SeqNado on a SLURM cluster with containerized environments:

```bash
seqnado pipeline atac --preset ss --queue short --scale-resources 1.0
```

The `ss` (SLURM Singularity) preset automatically configures job submission to SLURM with container isolation.

See CLI options for presets, queues, and scaling: [seqnado pipeline](cli.md#cli-seqnado-pipeline).

## SLURM Configuration

### Preset: `ss` (SLURM Singularity)

The SLURM Singularity preset uses:

- **Executor**: SLURM job scheduler
- **Containerization**: Singularity/Apptainer containers (recommended for HPC environments)
- **Default partition**: `short` (1 hour runtime, 3 GB memory per job)
- **Max parallel jobs**: 100

### Key Options

When running with the SLURM preset:

#### Queue/Partition
Specify the SLURM partition:

```bash
seqnado pipeline atac --preset ss --queue short
```

Common partitions:

- `short` - Default, 1 hour runtime, general purpose
- `long` - Extended runtime (varies by system)
- `gpu` - GPU-enabled partition (for GPU-accelerated steps like MCC)

#### Resource Scaling
Scale memory and runtime resources by a factor:

```bash
seqnado pipeline atac --preset ss --scale-resources 1.5
```

A factor of 1.5 increases all memory and time allocations by 50%. Use this if jobs are being killed due to resource limits.

#### Verbose Output
Print the generated Snakemake command before execution:

```bash
seqnado pipeline atac --preset ss --print-cmd
```

For a full list of flags, see [seqnado pipeline](cli.md#cli-seqnado-pipeline).

### GPU Support

For assays using GPU acceleration (e.g., MCC peak calling):

1. Submit to a GPU partition:
```bash
seqnado pipeline atac --preset ss --queue gpu
```

2. SeqNado automatically allocates GPU resources for GPU-enabled steps.

### Multiomics Workflows

For multiomics experiments combining multiple assays:

```bash
seqnado pipeline --preset ss --configfile config.yaml
```

Ensure your configuration file includes all assays and their respective genomes.

## Other Execution Options

While SLURM+Singularity is recommended for HPC, SeqNado provides several other profiles:

### Local Execution Presets

These options are useful for development, testing, or non-HPC environments:

#### `le` (Local Execution)
Standard local execution with environment isolation:

```bash
seqnado pipeline atac --preset le
```

Suitable for:

- Development machines
- Local testing
- Systems without job schedulers

#### `lc` (Local Cluster)
Local execution with job-level parallelization (GNU parallel):

```bash
seqnado pipeline atac --preset lc
```

Useful for:

- Multi-CPU workstations
- Quick local runs with parallelization

#### `ls` (Local Single-threaded)
Single-threaded local execution:

```bash
seqnado pipeline atac --preset ls
```

Useful for:

- Debugging
- Testing without parallelization

### Container Options

By default, all presets use Singularity/Apptainer containers for reproducibility. To use a local Conda environment instead:

Modify your workflow configuration YAML:

```yaml
container: null  # Disable containers
```

Then run with any preset. The workflow will execute with your current Python environment.

## Configuration Files

### Genome Configuration

Ensure genomes are configured before running pipelines:

```bash
seqnado genomes list atac
```

More details: [seqnado genomes](cli.md#cli-seqnado-genomes).

### Workflow Configuration

Generate an assay-specific configuration:

```bash
seqnado config atac --output config_atac.yaml
```

More details: [seqnado config](cli.md#cli-seqnado-config).

Edit the YAML to customize:

- Input/output paths
- Reference genomes
- Analysis parameters
- Resource allocations

### Design Metadata

Generate metadata CSV from FASTQ files:

```bash
seqnado design atac /path/to/fastqs/*.fastq.gz
```

More details: [seqnado design](cli.md#cli-seqnado-design).

## Troubleshooting

### Jobs Terminated Due to Resource Limits

If Snakemake jobs are killed with `out of memory` or timeout errors:

```bash
seqnado pipeline atac --preset ss --scale-resources 2.0 --queue long
```

- Increase `--scale-resources` factor
- Switch to a partition with longer walltime
- Check cluster documentation for memory limits

### Apptainer/Singularity Not Found

SeqNado prefers Apptainer (modern Singularity fork). If unavailable:

1. Install Apptainer on your cluster
2. Or switch to local environment execution (modify config YAML)
3. Or request Singularity installation from cluster admin

### SLURM Submission Errors

Verify SLURM is available:

```bash
sinfo  # List SLURM partitions
squeue # Check submitted jobs
```

Add `--print-cmd` to see the exact Snakemake command and submission details.

### Container Image Download Issues

SeqNado pulls container images on first use. If downloads fail due to network restrictions:

1. Use cluster's container cache if available
2. Pre-download and cache images on compute nodes
3. Contact cluster support for container registry access

## Advanced: Custom Snakemake Configuration

For advanced users, profiles can be customized by editing the configuration files in:

```
~/.config/seqnado/profiles/profile_slurm_singularity/
```

Refer to [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/) for configuration options.

## Support

For issues specific to cluster execution:

1. Check cluster documentation for scheduler syntax and available partitions
2. Verify containers are accessible from compute nodes
3. Test with `--preset le` on the login node first
4. Use `--print-cmd` to debug generated Snakemake commands
