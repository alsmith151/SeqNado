import ngs_pipeline.utils as utils

rule remove_blacklisted_regions:
    input: 
        bam = "aligned_and_filtered/{sample}.bam",
    output:
        sentinel = touch("flags/{sample}.blacklist.sentinel"),
    threads:
        1
    params:
        blacklist = config["blacklist"],
    log:
        "logs/blacklist/{sample}.log"
    run:

        if params.blacklist and os.path.exists(params.blacklist):

            cmd = f"""
                    bedtools intersect -v -b {params.blacklist} -a {input.bam} > {input.bam}.tmp &&
                    mv {input.bam}.tmp {input.bam} &&
                    echo "Removed blacklisted regions" > {log}
                   """

            if workflow.use_singularity:
                cmd = utils.get_singularity_command(command=cmd,
                                                    workflow=workflow,)
           
        else:
            cmd = f"""echo "No blacklisted regions specified" > {log}"""
        
        shell(cmd)


rule deeptools_shift_alignments:
    input:
        bam = "aligned_and_filtered/{sample}.bam",
    output:
        sentinel = touch("flags/{sample}.shifted.sentinel"),
    params:
        options = config["deeptools"]["alignmentsieve"]
    threads:
        8
    log:
        "logs/duplicate_removal/deeptools/{sample}.log",
    run:

        if config.get("shift_atac_reads"):
            
            cmd = """
            alignmentSieve -b {input.bam} -o {input.bam} -p {threads} {params.options} --ATACshift > {log} 2>&1 &&
            samtools sort -o {input.bam}.tmp {input.bam} -@ {threads} &&
            mv {input.bam}.tmp {input.bam} &&
            samtools index {input.bam} &&
            rm -f {input.bam}.tmp
            """

            if workflow.use_singularity:
                    cmd = utils.get_singularity_command(command=cmd,
                                                        workflow=workflow,)
            
        else:
            cmd = f"""echo "Will not shift reads" > {log}"""
        
        shell(cmd)


rule mark_filtering_complete:
    input:
        sentinel = "flags/{sample}.blacklist.sentinel",
        sentinel2 = "flags/{sample}.shifted.sentinel",
    output:
        sentinel = touch("flags/{sample}.filtering.complete.sentinel"),
    log:
        "logs/filtering/{sample}.log"
    shell:
        """echo "Filtering complete" > {log}"""
