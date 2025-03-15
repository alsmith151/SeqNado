from seqnado.helpers import check_options

rule methyldackel_extract:
    input:
        bam="seqnado_output/aligned/{sample}.bam",
    output:
        bdg = "seqnado_output/methylation/{sample}_CpG.bedGraph"
    params:
        fasta=config["fasta"],
        options=check_options(config["methyldackel"]["options"]),
        prefix="seqnado_output/methylation/{sample}"
    threads: config["methyldackel"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container:"library://cchahrou/seqnado/seqnado_meth.sif:latest"
    log:"seqnado_output/logs/methylation/methyldackel/{sample}.log"
    shell: """
        MethylDackel extract -@ {threads} {params.options} -o {params.prefix} {params.fasta} {input.bam} > {log} 2>&1
        """

rule methyldackel_bias:
    input:
        bam="seqnado_output/aligned/{sample}.bam"
    output:
        bias="seqnado_output/methylation/bias/{sample}.txt"
    params:
        fasta=config["fasta"],
        prefix="seqnado_output/methylation/bias/{sample}"
    threads: config["methyldackel"]["threads"]
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container:"library://cchahrou/seqnado/seqnado_meth.sif:latest"
    log:"seqnado_output/logs/methylation/methyldackel_bias/{sample}.log"
    shell: """
        MethylDackel mbias -@ {threads} --txt {params.fasta} {input.bam} {params.prefix} > {output.bias} 2> {log}
        """

rule taps_inverted:
    input:
        bdg="seqnado_output/methylation/{sample}_CpG.bedGraph"
    output:
        taps="seqnado_output/methylation/{sample}_CpG_inverted.bedGraph"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:"seqnado_output/logs/methylation/taps_inverted/{sample}.log"
    shell: """
        awk -v OFS="\t" '{{print $1, $2, $3, (100-$4), $5, $6}}' {input.bdg} > {output.taps} 2> {log}
        """