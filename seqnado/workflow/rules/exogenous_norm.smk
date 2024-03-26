from seqnado.design import NormGroups

NORM_GROUPS = NormGroups.from_design(DESIGN)

use rule align_paired as align_paired_spikein with:
    params:
        index=config["genome"]["indices"],
        options="--no-mixed --no-discordant",
    output:
        bam=temp("seqnado_output/aligned/spikein/raw/{sample}.bam"),
    resources:
        mem=lambda wildcards, attempt: f"{8 * 2 ** (attempt - 1)}GB",

use rule align_single as align_single_spikein with:
    output:
        bam=temp("seqnado_output/aligned/spikein/raw/{sample}.bam"),
    resources:
        mem=lambda wildcards, attempt: f"{2 * 2 ** (attempt - 1)}GB",


use rule sort_bam as sort_bam_spikein with:
    input:
        bam="seqnado_output/aligned/spikein/raw/{sample}.bam",
    output:
        bam=temp("seqnado_output/aligned/spikein/sorted/{sample}.bam"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}_sort.log",


use rule index_bam as index_bam_spikein with:
    input:
        bam=rules.sort_bam_spikein.output.bam,
    output:
        bai=temp("seqnado_output/aligned/spikein/sorted/{sample}.bam.bai"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}_index.log",


rule filter_bam_spikein:
    input:
        bam=rules.sort_bam_spikein.output.bam,
    output:
        bam=temp("seqnado_output/aligned/spikein/filtered/{sample}.bam"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}_filter.log",
    shell:
        """
    samtools view -b -F 3332 -q 30 -@ 8 {input.bam} > {output.bam} &&
    echo 'Filtered bam number of mapped reads:' > {log} 2>&1 &&
    samtools view -c {output.bam} >> {log} 2>&1
    """


use rule index_bam as index_bam_spikein_filtered with:
    input:
        bam=rules.filter_bam_spikein.output.bam,
    output:
        bai=temp("seqnado_output/aligned/spikein/filtered/{sample}.bam.bai"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}_filter_index.log",


rule split_bam:
    input:
        bam=rules.filter_bam_spikein.output.bam,
        bai=rules.index_bam_spikein_filtered.output.bai,
    output:
        ref_bam=temp("seqnado_output/aligned/spikein/{sample}_ref.bam"),
        exo_bam=temp("seqnado_output/aligned/spikein/{sample}_exo.bam"),
        stats="seqnado_output/aligned/spikein/{sample}_stats.tsv",
    params:
        genome_prefix=config["spikein_options"]["reference_genome"],
        exo_prefix=config["spikein_options"]["spikein_genome"],
        prefix="seqnado_output/aligned/spikein/{sample}",
        map_qual=30,
    log:
        "seqnado_output/logs/split_bam/{sample}.log",
    shell:
        """
        samtools view -h {input.bam} | awk '{{if($0 ~ /^@/ || $3 ~ /^chr/) print}}' | samtools view -b -o {output.ref_bam} &&
        samtools view -h {input.bam} | awk '{{if($0 ~ /^@/ || $3 ~ /^{params.exo_prefix}/) print}}' | samtools view -b -o {output.exo_bam} &&
        echo -e "sample\treference_reads\tspikein_reads" > {output.stats} &&
        echo -e "{wildcards.sample}\t$(samtools view -c {output.ref_bam})\t$(samtools view -c {output.exo_bam})" >> {output.stats}
        """


rule move_ref_bam:
    input:
        bam=rules.split_bam.output.ref_bam,
    output:
        bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
    shell:
        """
    mv {input.bam} {output.bam}
    """


if config["spikein_options"]["normalisation_method"] == "orlando":

    rule calculate_normalisation_factors:
        input:
            lambda wc: expand(rules.split_bam.output.stats, sample=[sample for sample in NORM_GROUPS.get_grouped_samples(wc.group)]),
        output:
            normalisation_table="seqnado_output/resources/{group}_normalisation_factors.tsv",
            normalisation_factors="seqnado_output/resources/{group}_normalisation_factors.json",
        log:
            "seqnado_output/logs/normalisation_factors_{group}.log",
        script:
            "../scripts/calculate_spikein_norm_orlando.py"

elif config["spikein_options"]["normalisation_method"] == "with_input":

    rule calculate_normalisation_factors:
        input:
            lambda wc: expand(rules.split_bam.output.stats, sample=[sample for sample in NORM_GROUPS.get_grouped_samples(wc.group)]),
        output:
            normalisation_table="seqnado_output/resources/{group}_normalisation_factors.tsv",
            normalisation_factors="seqnado_output/resources/{group}_normalisation_factors.json",
        log:
            "seqnado_output/logs/normalisation_factors_{group}.log",
        script:
            "../scripts/calculate_spikein_norm_factors.py"


ruleorder: align_paired_spikein > align_single_spikein
