[← Back to main page](index.md)

# Outputs

All SeqNado analysis results are organized within the `seqnado_output/` directory (or your custom output directory). This page describes the general structure and types of files you can expect.

## General Structure

The output directory follows a structured organization:

- `qc/`: Quality control reports (FastQC, MultiQC, etc.)
- `alignments/`: BAM files and indexing information.
- `peaks/`: (Assay-specific) Peak calling results (BED, BEDGraph, peaks files).
- `tracks/`: BigWig tracks for visualization in genome browsers.
- `plots/`: Diagnostic plots, heatmaps, and TSS enrichment plots.
- `report/`: Summary reports of the entire pipeline run.
- `geo/`: (Optional) Prepared metadata and symlinks for GEO submission.

## Assay-Specific Outputs

### ATAC-seq and ChIP-seq

- **Peaks**: Found in `peaks/` using MACS2 or other specified callers.
- **Coverage**: BigWig tracks generated in `tracks/`, often scaled or normalized.
- **QC**: TSS enrichment scores and fragment size distributions in `qc/`.

### RNA-seq

- **Quantification**: Gene and transcript level counts in `quant/`.
- **Alignments**: Genomic coordinates BAM files in `alignments/`.

### CUT&Tag (CAT)

- **Peaks**: Peaks in `peaks/` with SEACR or configured caller.
- **Coverage**: BigWigs in `tracks/` optimized for small fragments.
- **QC**: Low-input focused metrics (e.g., fragment size summaries) in `qc/`.

### SNP Analysis (SNP)

- **Variants**: VCF/BCF files and summaries in `variants/` (if enabled).
- **Alignments**: Sorted/indexed BAMs in `alignments/`.
- **QC**: Mapping/variant-calling metrics in `qc/`.

### Methylation (METH)

- **Methylation Calls**: Cytosine/CpG reports (e.g., BED/bedGraph) in `methylation/`.
- **Alignments**: BAMs in `alignments/` (bisulfite-aware processing).
- **QC**: Conversion rates and summary metrics in `qc/`.

### MCC (MCC)

- **Contacts**: Contact maps/matrices (e.g., cooler or equivalent formats) in `interactions/`.
- **Peaks/Enrichment**: MCC-specific enrichment or peak outputs in `peaks/` (if configured).
- **QC**: MCC run metrics and summaries in `qc/`.

### CRISPR Screens (CRISPR)

- **Guide Counts**: Count tables per guide/target in `quant/` or `counts/`.
- **Screen Analysis**: Differential hit results (e.g., MAGeCK outputs) in `report/` or `analysis/`.
- **QC**: Library representation and mapping metrics in `qc/`.

### Multiomics (MULTIOMICS)

- **Per‑Assay Outputs**: All of the above, organized under a common project with per‑assay subfolders.
- **Integrated Reports**: Cross‑assay summaries in `report/` and combined plots in `plots/`.

## Summary Report

The primary entry point for exploring your results is the `seqnado_report.html` found in the `report/` folder. It provides an interactive overview of all samples and QC metrics.

To rerun or adjust resources/presets, see the CLI reference: [seqnado pipeline](cli.md#cli-seqnado-pipeline).
