checkpoint methylation_split_bams:
    input:
        OUTPUT_DIR + "/aligned/{sample}.bam"
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_{genome}.bam"),
        bai=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_{genome}.bam.bai")
    params:
        ref_genome=config["genome"]["name"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/methylation/split_bams/{sample}_{genome}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation_split_bams/{sample}_{genome}.tsv"
    message: "Splitting BAM for sample {wildcards.sample} into genome {wildcards.genome}"
    shell: """
    if [[ "{wildcards.genome}" == "{params.ref_genome}" ]]; then
        samtools view -h {input} | awk '{{if($0 ~ /^@/ || $3 ~ /^chr/) print}}' | samtools view -b -o {output.bam} 2> {log}
    else
        samtools view -h {input} | awk '{{if($0 ~ /^@/ || $3 ~ /^{wildcards.genome}/) print}}' | samtools view -b -o {output.bam} 2> {log}
    fi
    samtools index {output.bam}
    """


def get_split_bam(wildcards):
    """Retrieve the resolved BAM file from the checkpoint."""
    checkpoint_output = checkpoints.methylation_split_bams.get(sample=wildcards.sample, genome=wildcards.genome).output
    return checkpoint_output[0] 


rule methyldackel_bias:
    input:
        bam=get_split_bam
    output:
        bias=OUTPUT_DIR + "/methylation/methyldackel/bias/{sample}_{genome}.txt",
    params:
        fasta=CONFIG.genome.fasta,
        prefix=OUTPUT_DIR + "/methylation/methyldackel/bias/{sample}_{genome}"
    threads: CONFIG.third_party_tools.methyldackel.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container:"oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:OUTPUT_DIR + "/logs/methylation/methyldackel/bias/{sample}_{genome}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation/methyldackel/bias/{sample}_{genome}.tsv"
    message: "Running MethylDackel bias for sample {wildcards.sample} and genome {wildcards.genome}"
    shell: """
    MethylDackel mbias -@ {threads} --txt {params.fasta} {input.bam} {params.prefix} > {output.bias} 2> {log}
    """

rule calculate_conversion:
    input:
        bias=expand(OUTPUT_DIR + "/methylation/methyldackel/bias/{sample}_{genome}.txt", sample=SAMPLE_NAMES, genome=SPIKEIN_GENOMES)
    output:
        conversion=OUTPUT_DIR + "/methylation/methylation_conversion.tsv",
        plot=OUTPUT_DIR + "/methylation/methylation_conversion.png"
    params:
        assay=CONFIG.methylation.method,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/methylation/conversion.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation/calculate_conversion.tsv"
    message: "Calculating methylation conversion rates across all samples"
    script: "../scripts/methylation_conversion.py"
    

rule methyldackel_extract:
    input:
        bam=get_split_bam
    output:
        bdg=OUTPUT_DIR + "/methylation/methyldackel/{sample}_{genome}_CpG.bedGraph"
    params:
        fasta=CONFIG.genome.fasta,
        options=CONFIG.third_party_tools.methyldackel.options,
        prefix=OUTPUT_DIR + "/methylation/methyldackel/{sample}_{genome}"
    threads: CONFIG.third_party_tools.methyldackel.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/methylation/methyldackel/{sample}_{genome}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation/methyldackel/{sample}_{genome}.tsv"
    message: "Running MethylDackel extract for sample {wildcards.sample} and genome {wildcards.genome}"
    shell: """
    MethylDackel extract -@ {threads} {params.options} -o {params.prefix} {params.fasta} {input.bam} > {log} 2>&1
    """


rule taps_inverted:
    input:
        bdg=OUTPUT_DIR + "/methylation/methyldackel/{sample}_{genome}_CpG.bedGraph"
    output:
        taps=OUTPUT_DIR + "/methylation/methyldackel/{sample}_{genome}_CpG_TAPS.bedGraph"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/methylation/taps_inverted/{sample}_{genome}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation/taps_inverted/{sample}_{genome}.tsv"
    message: "Converting to TAPS methylation for sample {wildcards.sample} and genome {wildcards.genome}"
    shell: """
    awk -v OFS="\t" '{{print $1, $2, $3, (100-$4), $5, $6}}' {input.bdg} > {output.taps} 2> {log}
    rm {input.bdg}
    """

localrules:
    calculate_conversion
