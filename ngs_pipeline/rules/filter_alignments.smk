import ngs_pipeline.utils as utils

rule remove_blacklisted_regions:
    input: 
        bam = "aligned_and_filtered/{sample}.bam",
    output:
        log = "logs/blacklist/{sample}.log"
    threads:
        1
    params:
        blacklist = config["blacklist"],
    log:
        "logs/blacklisted_regions_removed_{sample}.log"
    run:

        if params.blacklist and os.path.exists(params.blacklist):

            cmd = f"""
                    bedtools intersect -v -b {params.blacklist} -a {input.bam} > {input.bam}.tmp &&
                    mv {input.bam}.tmp {input.bam} &&
                    echo "Removed blacklisted regions" > {output.log}
                   """

            if workflow.use_singularity:
                cmd = utils.get_singularity_command(command=cmd,
                                                    workflow=workflow,)
           
        else:
            cmd = f"""echo "No blacklisted regions specified" > {output.log}"""
        
        shell(cmd)
