from seqnado.helpers import check_options

checkpoint methylation_split_bams:
    input:
        "seqnado_output/aligned/{sample}.bam"
    output:
        bam=temp("seqnado_output/aligned/spikein/{sample}_{genome}.bam"),
        bai=temp("seqnado_output/aligned/spikein/{sample}_{genome}.bam.bai")
    params:
        ref_genome=config["genome"]["name"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/methylation/split_bams/{sample}_{genome}.log"
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
        bias="seqnado_output/methylation/methyldackel/bias/{sample}_{genome}.txt",
        OB_svg="seqnado_output/methylation/methyldackel/bias/{sample}_{genome}_OB.svg",
        OT_svg="seqnado_output/methylation/methyldackel/bias/{sample}_{genome}_OT.svg"
    params:
        fasta=config["fasta"],
        prefix="seqnado_output/methylation/methyldackel/bias/{sample}_{genome}"
    threads: config["methyldackel"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container:"library://cchahrou/seqnado/seqnado_meth.sif:latest"
    log:"seqnado_output/logs/methylation/methyldackel/bias/{sample}_{genome}.log"
    shell: """
        MethylDackel mbias -@ {threads} --txt {params.fasta} {input.bam} {params.prefix} > {output.bias} 2> {log}
    """

rule calculate_conversion:
    input:
        bias=expand("seqnado_output/methylation/methyldackel/bias/{sample}_{genome}.txt", sample=SAMPLE_NAMES, genome=SPIKEIN_GENOMES)
    output:
        conversion="seqnado_output/methylation/methylation_conversion.tsv",
        plot="seqnado_output/methylation/methylation_conversion.png"
    params:
        assay=config["methylation_assay"],
    log:
        "seqnado_output/logs/methylation/conversion.log"
    script: "../scripts/methylation_conversion.py"
    

rule methyldackel_extract:
    input:
        bam=get_split_bam
    output:
        bdg="seqnado_output/methylation/methyldackel/{sample}_{genome}_CpG.bedGraph"
    params:
        fasta=config["fasta"],
        options=check_options(config["methyldackel"]["options"]),
        prefix="seqnado_output/methylation/methyldackel/{sample}_{genome}"
    threads: config["methyldackel"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "library://cchahrou/seqnado/seqnado_meth.sif:latest"
    log:
        "seqnado_output/logs/methylation/methyldackel/{sample}_{genome}.log"
    shell: """
        MethylDackel extract -@ {threads} {params.options} -o {params.prefix} {params.fasta} {input.bam} > {log} 2>&1
    """


rule taps_inverted:
    input:
        bdg="seqnado_output/methylation/methyldackel/{sample}_{genome}_CpG.bedGraph"
    output:
        taps="seqnado_output/methylation/methyldackel/{sample}_{genome}_CpG_TAPS.bedGraph"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/methylation/taps_inverted/{sample}_{genome}.log"
    shell: """
        awk -v OFS="\t" '{{print $1, $2, $3, (100-$4), $5, $6}}' {input.bdg} > {output.taps} 2> {log}
        rm {input.bdg}
    """

localrules:
    calculate_conversion
