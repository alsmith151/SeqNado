rule ignore_duplicates:
    input:
        bam="aligned/{sample}.bam",
        index = "aligned/{sample}.bam.bai",
    output:
        bam= "aligned_and_filtered/{sample}.bam",
        log = "logs/duplicate_removal/not_removed/{sample}.log",
    log:
        "logs/duplicate_removal/not_removed/{sample}.log",
    run:
        abspath = os.path.abspath(input.bam)
        shell(f"""ln -s {abspath} {output.bam} && 
                    ln -s {input.bam}.bai {output.bam}.bai
                """)