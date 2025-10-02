from seqnado.helpers import define_memory_requested, define_time_requested
from seqnado.config.third_party_tools import CommandLineArguments
from seqnado import Assay
import pandas as pd
import re


def format_feature_counts(counts: str) -> pd.DataFrame:
    counts = pd.read_csv(counts, sep="\t", comment="#")
    counts = counts.set_index("Geneid")
    counts = counts.drop(
        columns=["Chr", "Start", "End", "Strand", "Length"], errors="ignore"
    )
    count_files = get_count_files(None)
    counts = counts.loc[:, count_files]
    return counts


def create_metadata(counts: pd.DataFrame) -> pd.DataFrame:
    count_files = get_count_files(None)
    counts = counts.loc[:, count_files]
    df = counts.columns.to_frame("sample_name")
    return df


def get_scaling_factor(wildcards, scale_path: str) -> float:
    df = pd.read_csv(scale_path, sep="\t", index_col=0)
    try:
        factor = df.loc[wildcards.sample, "norm.factors"]
    except KeyError:
        factor = 1

    return factor


def _sample_group(mapping_name: str, wildcards) -> str:
    groups = SAMPLE_GROUPINGS.get_grouping(mapping_name)
    return groups.sample_to_group()[wildcards.sample]


def get_norm_factor_spikein(wildcards, negative=False):
    import json

    group = _sample_group("normalisation", wildcards)
    with open(f"seqnado_output/resources/{group}_normalisation_factors.json") as f:
        norm_factors = json.load(f)

    if not negative:
        return norm_factors[wildcards.sample]
    else:
        return -norm_factors[wildcards.sample]


def format_deeptools_bamcoverage_options(wildcards):
    options_str = str(
        CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments
    )

    # Remove normalization flags; we manage scale explicitly
    options_str = re.sub(r"--normalizeUsing [a-zA-Z]+", "", options_str)
    options_str = re.sub(r"--scaleFactor [0-9.]+", "", options_str)

    # Remove extendReads flags for single-end
    if not INPUT_FILES.is_paired_end(wildcards.sample):
        options_str = re.sub(r"--extendReads", "", options_str)
        options_str = re.sub(r"\s-e(\s|$)", " ", options_str)

    return options_str.strip()


def format_homer_make_bigwigs_options(wildcards):
    norm = int(get_norm_factor_spikein(wildcards) * 1e7)
    options_str = str(
        CONFIG.third_party_tools.homer.make_bigwig.command_line_arguments
    )

    # Ensure any existing -norm is removed, then add ours
    options_str = re.sub(r"-norm [0-9.]+", "", options_str)
    options_str = f"{options_str} -norm {norm}".strip()
    return options_str


# CSAW Method
rule tile_regions:
    input:
        chromsizes=CONFIG.genome.chromosome_sizes,
    output:
        genome_tiled="seqnado_output/resources/genome_tiled.gtf",
    params:
        tile_size=lambda wc: CONFIG.genome.bin_size if CONFIG.genome.bin_size is not None else 10000,
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    script:
        "../../scripts/tile_genome.py"


def get_count_files(wildcards):
    import pathlib

    files = []

    if ASSAY in [Assay.CHIP, Assay.CAT]:
        df = INPUT_FILES.to_dataframe()
        for row in df.itertuples():
            try:
                f1 = pathlib.Path(row.ip_r1)
                f2 = pathlib.Path(row.ip_r2)
                if f1.exists() and f2.exists():
                    files.append(
                        f"seqnado_output/aligned/{row.sample_name}_{row.ip}.bam"
                    )
            except Exception:
                pass

    elif ASSAY == Assay.ATAC:
        df = INPUT_FILES.to_dataframe()
        files = expand("seqnado_output/aligned/{sample}.bam", sample=df.sample_name)

    return files


rule count_bam:
    input:
        bam=lambda wc: get_count_files(wc),
        tiles="seqnado_output/resources/genome_tiled.gtf",
    output:
        counts="seqnado_output/counts/counts.tsv",
    params:
        options="-p --countReadPairs",
    log:
        "seqnado_output/logs/counts/readcounts.log",
    threads: 8
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        "featureCounts -a {input.tiles} -a {input.tiles} -t tile -o {output.counts} {input.bam} -T {threads} {params.options} > {log} 2>&1"


rule setup_for_scaling_factors:
    input:
        counts="seqnado_output/counts/counts.tsv",
    output:
        formatted_counts="seqnado_output/counts/{group}_formatted_counts.tsv",
        metadata="seqnado_output/counts/{group}_metadata.tsv",
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    run:
        counts = format_feature_counts(input[0])
        counts.to_csv(output[0], sep="\t", index=False)

        metadata = create_metadata(counts)
        metadata.to_csv(output[1], sep="\t", index=False, header=False)


rule calculate_scaling_factors:
    input:
        formatted_counts="seqnado_output/counts/{group}_formatted_counts.tsv",
        metadata="seqnado_output/counts/{group}_metadata.tsv",
    output:
        scaling_factors="seqnado_output/resources/{group}_scaling_factors.tsv",
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    script:
        "../../scripts/calculate_scaling_factors.R"


rule calculate_scaling_factors_spikein:
    input:
        counts="seqnado_output/readcounts/feature_counts/read_counts.tsv",
        metadata="seqnado_output/design.csv",
    output:
        size_factors="seqnado_output/resources/all_normalisation_factors.json"
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    params:
        spikein_genes=["AmpR_seq", "Cas9_5p_seq", "Cas9_3p_seq"],
    log:
        "seqnado_output/logs/normalisation_factors.log"
    script:
        "../../scripts/calculate_spikein_norm_factors_rna.R"


rule deeptools_make_bigwigs_scale:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: f"seqnado_output/resources/{_sample_group('scaling', wc)}_scaling_factors.tsv",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/csaw/{sample}.bigWig",
    params:
        scale=lambda wc: get_scaling_factor(
            wc,
            f"seqnado_output/resources/{_sample_group('scaling', wc)}_scaling_factors.tsv",
        ),
        options=lambda wc: format_deeptools_bamcoverage_options(wc)
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),    
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/pileups/deeptools/scaled/{sample}.log",
    shell:
        "bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options} > {log} 2>&1"


use rule deeptools_make_bigwigs_scale as deeptools_make_bigwigs_spikein with:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: f"seqnado_output/resources/{_sample_group('normalisation', wc)}_normalisation_factors.json",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/spikein/{sample}.bigWig",
    params:
        options=lambda wildcards: format_deeptools_bamcoverage_options(wildcards),
        scale=get_norm_factor_spikein,


rule deeptools_make_bigwigs_rna_spikein_plus:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: f"seqnado_output/resources/{_sample_group('normalisation', wc)}_normalisation_factors.json",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/spikein/{sample}_plus.bigWig",
    params:
        options=lambda wildcards: format_deeptools_bamcoverage_options(wildcards),
        scale=get_norm_factor_spikein,
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/pileups/deeptools/spikein/{sample}_plus.log",
    shell:
        "bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} --scaleFactor {params.scale} {params.options} --filterRNAstrand forward > {log} 2>&1"


rule deeptools_make_bigwigs_rna_spikein_minus:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
        bai="seqnado_output/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: f"seqnado_output/resources/{_sample_group('normalisation', wc)}_normalisation_factors.json",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/spikein/{sample}_minus.bigWig",
    params:
        options=lambda wildcards: format_deeptools_bamcoverage_options(wildcards),
        scale=lambda wc: get_norm_factor_spikein(wc, negative=True),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/pileups/deeptools/spikein/{sample}_minus.log",
    shell:
        "bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} --scaleFactor {params.scale} {params.options} --filterRNAstrand reverse > {log} 2>&1"
