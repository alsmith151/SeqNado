
if REMOVE_DUPLICATES_DEEPTOOLS or FILTER_ALIGNMENTS:
    rule deeptools_filter_alignments:
        input:
            bam = "aligned/{sample}.bam",
        output:
            bam = "aligned_and_filtered/{sample}.bam",
            index = "aligned_and_filtered/{sample}.bam.bai",
            log = "logs/deeptools_filter_alignments/{sample}.log",
        params:
            threads = config["pipeline"]["n_cores"],
            deduplicate = "--ignoreDuplicates" if REMOVE_DUPLICATES else "",
            options = config["alignments"]["filter_options"] if config["alignments"]["filter_options"] else "",
        log:
            "logs/deeptools_filter_alignments_{sample}.log"
        shell:
            """
            alignmentSieve -b {input.bam} -o {output.bam} -p {params.threads} {params.deduplicate} {params.options} 2>&1 &&
            samtools sort -o {output.bam}.tmp {output.bam} -@ {params.threads} &&
            mv {output.bam}.tmp {output.bam} &&
            samtools index {output.bam} &&
            rm -f {output.bam}.tmp
            """
if REMOVE_DUPLICATES_PICARD:
    rule mark_duplicates:
        input:
            bams="aligned/{sample}.bam",
        # optional to specify a list of BAMs; this has the same effect
        # of marking duplicates on separate read groups for a sample
        # and then merging
        output:
            bam="aligned_and_filtered/{sample}.bam",
            metrics="aligned_and_filtered/{sample}.metrics.txt",
            log = "logs/picard_mark_duplicates_{sample}.log"
        log:
            "logs/picard_mark_duplicates_{sample}.log",
        params:
            extra="--REMOVE_DUPLICATES true",
        # optional specification of memory usage of the JVM that snakemake will respect with global
        # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
        # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
        # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
        resources:
            mem_mb=1024 * 4,
        wrapper:
            "v1.5.0/bio/picard/markduplicates"
        
    rule index_bam:
        input:
            bam="aligned_and_filtered/{sample}.bam",
        output:
            index="aligned_and_filtered/{sample}.bam.bai",
        params:
            threads=config["pipeline"]["n_cores"],
        shell:
            "samtools index {input.bam} -@ {params.threads}"

if not REMOVE_DUPLICATES:

    rule ignore_duplicates:
        input:
            bam="aligned/{sample}.bam",
            index = "aligned/{sample}.bam.bai",
        output:
            bam= "aligned_and_filtered/{sample}.bam",
            log = "logs/duplicates_ignored_{sample}.log",
        log:
            "logs/duplicates_ignored_{sample}.log",
        run:
            abspath = os.path.abspath(input.bam)
            shell(f"""ln -s {abspath} {output.bam} && 
                    ln -s {input.bam}.bai {output.bam}.bai""")



rule remove_blacklisted_regions:
    input: 
        bam = "aligned_and_filtered/{sample}.bam",
    output:
        log = "logs/blacklisted_regions_removed_{sample}.log"
    params:
        blacklist = config["genome"]["blacklist"],
        threads = config["pipeline"]["n_cores"],
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
