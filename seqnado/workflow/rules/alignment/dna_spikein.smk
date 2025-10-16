
###################################################
# Spike-in alignment and processing for DNA assays #
###################################################

use rule align_paired as align_paired_spikein with:
    params:
        options="--no-mixed --no-discordant",
        index=CONFIG.genome.index.prefix,
        rg="--rg-id {sample} --rg SM:{sample}",
    output:
        bam=temp("seqnado_output/aligned/spikein/raw/{sample}.bam"),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"


use rule align_single as align_single_spikein with:
    output:
        bam=temp("seqnado_output/aligned/spikein/raw/{sample}.bam"),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"


use rule sort_bam as sort_bam_spikein with:
    input:
        bam="seqnado_output/aligned/spikein/raw/{sample}.bam",
    output:
        bam=temp("seqnado_output/aligned/spikein/sorted/{sample}.bam"),
        read_log=temp("seqnado_output/aligned/spikein/sorted/{sample}_read.log"),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/aligned_spikein/{sample}_sort.log",


use rule index_bam as index_bam_spikein with:
    input:
        bam=rules.sort_bam_spikein.output.bam,
    output:
        bai=temp("seqnado_output/aligned/spikein/sorted/{sample}.bam.bai"),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/aligned_spikein/{sample}_index.log",


rule filter_bam_spikein:
    input:
        bam=rules.sort_bam_spikein.output.bam,
    output:
        bam=temp("seqnado_output/aligned/spikein/filtered/{sample}.bam"),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        "seqnado_output/logs/aligned_spikein/{sample}_filter.log",
    shell:
        """
    samtools view -b -F 260 -@ 8 {input.bam} > {output.bam} &&
    echo 'Filtered bam number of mapped reads:' > {log} 2>&1 &&
    samtools view -c {output.bam} >> {log} 2>&1
    """


use rule index_bam as index_bam_spikein_filtered with:
    input:
        bam=rules.filter_bam_spikein.output.bam,
    output:
        bai=temp("seqnado_output/aligned/spikein/filtered/{sample}.bam.bai"),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
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
        genome_prefix=config.get("spikein_options", {}).get("reference_genome"),
        exo_prefix=config.get("spikein_options", {}).get("spikein_genome"),
        prefix="seqnado_output/aligned/spikein/{sample}",
        map_qual=30,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
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
    container: None
    shell:
        """
    mv {input.bam} {output.bam}
    """

