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
            
            shell("""
            bedtools intersect -v -b {params.blacklist} -a {input.bam} > {input.bam}.tmp &&
            mv {input.bam}.tmp {input.bam} &&
            echo "Removed blacklisted regions" > {output.log}
            """)
        else:
            shell("""echo "No blacklisted regions specified" > {output.log}""")
