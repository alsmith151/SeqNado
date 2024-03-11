from seqnado import utils

def format_feature_counts(counts: str) -> pd.DataFrame:
    counts = pd.read_csv(input.counts, sep="\t", comment="#")
    counts = counts.set_index("Geneid")
    counts = counts.drop(columns=["Chr", "Start", "End", "Strand", "Length"], errors="ignore")
    return counts

def create_metadata(counts: pd.DataFrame) -> pd.DataFrame:
    return counts.columns.str.replace(".bam", "")

def get_scaling_factor(wildcards: Any,  scale_path: str) -> float:
    df = pd.read_csv(scale_path, sep="\t", header=None, index_col=0)
    return df.loc[wildcards.sample, "norm.factors"]


def get_norm_factor_spikein(wildcards):
    with open(f"seqnado_output/resources/{wildcards.group}_spikein_factors.json") as f:
        norm_factors = json.load(f)
    return norm_factors[wildcards.sample]


def format_deeptools_bamcoverage_options(wildcards, invert_scale: bool = False):
    import re

    norm = get_norm_factor_spikein(wildcards)
    norm = norm if invert_scale else -norm


    options = utils.check_options(config["deeptools"]["bamcoverage"])

    if "--normalizeUsing" in options:
        options = re.sub("--normalizeUsing [a-zA-Z]+", "", options)

    if "--scaleFactor" in options:
        options = re.sub("--scaleFactor [0-9.]+", "", options)

    options += f" --scaleFactor {norm}"

    return options


def format_homer_make_bigwigs_options(wildcards):
    import re

    norm = int(get_norm_factor_spikein(wildcards) * 1e7)

    options = utils.check_options(config["homer"]["makebigwig"])

    if "-norm" in options:
        options = re.sub("-scale [0-9.]+", "", options)

    options += f" -norm {norm}"

    return options

# CSAW Method
rule tile_regions:
    input:
        chromsizes=config["chromsizes"],
    output:
        genome_tiled="seqnado_output/resources/genome_tiled.gtf"
    params:
        tile_size=config["tile_size"]
    run:
        import pyranges as pr
        chromsizes = (
            pd.read_csv(chromsizes, sep="\t", header=None).set_index(0)[1].to_dict()
        )
        genome_tiled = pr.gf.tile_genome(chromsizes, tile_size=tile_size)
        genome_tiled = genome_tiled.df.assign(feature="tile", gene_id=lambda df: df.index.astype(str)).pipe(pr.PyRanges)
        genome_tiled.to_gtf(output.genome_tiled)


rule count_bam:
    input:
        bam=expand("seqnado_output/aligned/{sample}.bam", sample=DESIGN.sample_names),
        tiles="seqnado_output/resources/genome_tiled.gtf"
    output:
        counts="seqnado_output/counts/counts.tsv"
    threads: 8
    shell:
        "featureCounts -a {input.tiles} -a {input.tiles} -t tiles -o {output.counts} {input.bam} -T {threads} -p --countReadPairs"

rule setup_for_scaling_factors:
    input:
        counts="seqnado_output/counts/{group}_counts.tsv"
    output:
        formatted_counts="seqnado_output/counts/{group}_formatted_counts.tsv",
        metadata="seqnado_output/counts/{group}_metadata.tsv",
    run:
        counts = format_counts(input.counts)
        counts.to_csv(output.formatted_counts, sep="\t")
        
        metadata = create_metadata(counts)
        metadata.to_csv(output.metadata, sep="\t", index=False, header=False)


rule calculate_scaling_factors:
    input:
        formatted_counts="seqnado_output/counts/{group}_formatted_counts.tsv",
        metadata="seqnado_output/counts/{group}_metadata.tsv",
    output:
        scaling_factors="seqnado_output/resources/{group}_scaling_factors.tsv"
    script:
        "../scripts/calculate_scaling_factors.R"
           

rule deeptools_make_bigwigs_scale:
    input:
        bam="seqnado_output/bams/{sample}.bam",
        bai="seqnado_output/bams/{sample}.bam.bai",
        scaling_factors=lambda wc: f"seqnado_output/resources/{utils.get_group_for_sample(wc)}_scaling_factors.tsv",
    output:
        bigwig="seqnado_output/bigwigs/window-norm/{sample}.bigWig",
    params:
        scale=lambda wc: get_scaling_factor(wc, f"seqnado_output/resources/{utils.get_group_for_sample(wc)}_scaling_factors.tsv"),
        options=utils.check_options(config["deeptools"]["bamcoverage"]),
    threads: 8
    shell:
        "bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options}"
        


use rule deeptools_make_bigwigs_scale as deeptools_make_bigwigs_spikein with:
    input:
        scaling_factors=lambda wc: f"seqnado_output/resources/{utils.get_group_for_sample(wc)}_spikein_factors.json",
    output:
        bigwig="seqnado_output/bigwigs/spikein-norm/{sample}.bigWig",
    params:
        options=lambda wildcards: format_deeptools_bamcoverage_options(wildcards),


rule deeptools_make_bigwigs_rna_spikein_plus:
    input:
        bam="seqnado_output/bams/{sample}.bam",
        bai="seqnado_output/bams/{sample}.bam.bai",
        scaling_factors=lambda wc: f"seqnado_output/resources/{utils.get_group_for_sample(wc)}_spikein_factors.json",
    output:
        bigwig="seqnado_output/bigwigs/spikein-norm/{sample}_plus.bigWig",
    params:
        options=lambda wildcards: format_deeptools_bamcoverage_options(wildcards),
    threads: 8
    shell:
        "bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} {params.options} --filterRNAstrand forward"

rule deeptools_make_bigwigs_rna_spikein_minus:
    input:
        bam="seqnado_output/bams/{sample}.bam",
        bai="seqnado_output/bams/{sample}.bam.bai",
        scaling_factors=lambda wc: f"seqnado_output/resources/{utils.get_group_for_sample(wc)}_spikein_factors.json",
    output:
        bigwig="seqnado_output/bigwigs/spikein-norm/{sample}_minus.bigWig",
    params:
        options=lambda wildcards: format_deeptools_bamcoverage_options(wildcards, invert_scale=True),
    threads: 8
    shell:
        "bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} {params.options} --filterRNAstrand reverse"
    