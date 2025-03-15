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

