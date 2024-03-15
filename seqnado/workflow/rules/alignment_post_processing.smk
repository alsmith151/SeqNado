from seqnado.helpers import check_options


rule sort_bam:
    input:
        bam="seqnado_output/aligned/raw/{sample}.bam",
    output:
        bam=temp("seqnado_output/aligned/sorted/{sample}.bam"),
    resources:
        mem=lambda wildcards, attempt: 4000 * 2**attempt,
    threads: 8
    log:
        "seqnado_output/logs/sorted/{sample}.log",
    shell:
        """
        samtools sort {input.bam} -@ {threads} -o {output.bam} -m 900M &&
        echo 'Sorted bam number of mapped reads:' > {log} 2>&1 &&
        samtools view -f 2 -c {output.bam} >> {log} 2>&1
        """


rule index_bam:
    input:
        bam="seqnado_output/aligned/sorted/{sample}.bam",
    output:
        bai=temp("seqnado_output/aligned/sorted/{sample}.bam.bai"),
    threads: 1
    resources:
        mem=1000,
    shell:
        "samtools index -@ {threads} -b {input.bam}"


if config["remove_blacklist"] and os.path.exists(config.get("blacklist", "")):

    rule remove_blacklisted_regions:
        input:
            bam="seqnado_output/aligned/sorted/{sample}.bam",
            bai=rules.index_bam.output.bai,
        output:
            bam=temp("seqnado_output/aligned/blacklist_regions_removed/{sample}.bam"),
            bai=temp(
                "seqnado_output/aligned/blacklist_regions_removed/{sample}.bam.bai"
            ),
        threads: 1
        params:
            blacklistcheck_options(config["blacklist"]),
        resources:
            mem="5GB",
            runtime="4h",
        log:
            "seqnado_output/logs/blacklist/{sample}.log",
        shell:
            """
            bedtools intersect -v -b {params.blacklist} -a {input.bam} > {output.bam} &&
            samtools index -b {output.bam} -o {output.bai} &&
            echo "Removed blacklisted regions" > {log} &&
            echo 'Number of mapped reads' >> {log} 2>&1 &&
            samtools view -f 2 -c {output.bam} >> {log} 2>&1
            """

else:

    rule ignore_blacklisted_regions:
        input:
            bam="seqnado_output/aligned/sorted/{sample}.bam",
            bai=rules.index_bam.output.bai,
        output:
            bam=temp("seqnado_output/aligned/blacklist_regions_removed/{sample}.bam"),
            bai=temp(
                "seqnado_output/aligned/blacklist_regions_removed/{sample}.bam.bai"
            ),
        threads: 1
        resources:
            mem="1GB",
        log:
            "seqnado_output/logs/blacklist/{sample}.log",
        shell:
            """
                mv {input.bam} {output.bam} &&
                mv {input.bai} {output.bai} &&
                echo "No blacklisted regions specified" > {log} &&
                echo 'Number of mapped reads' >> {log} 2>&1 &&
                samtools view -f 2 -c {output.bam} >> {log} 2>&1
                """


if config["remove_pcr_duplicates_method"] == "picard":

    rule remove_duplicates_using_picard:
        input:
            bam="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam",
            bai="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam.bai",
        output:
            bam=temp("seqnado_output/aligned/duplicates_removed/{sample}.bam"),
            bai=temp("seqnado_output/aligned/duplicates_removed/{sample}.bam.bai"),
            metrics=temp("seqnado_output/aligned/duplicates_removed/{sample}.metrics"),
        threads: 8
        params:
            optionscheck_options(config["picard"]["options"]),
        resources:
            mem="5GB",
            runtime="4h",
        log:
            "seqnado_output/logs/duplicates/{sample}.log",
        shell:
            """
            picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} --REMOVE_DUPLICATES true --CREATE_INDEX true {params.options} > {log} 2>&1 &&
            mv seqnado_output/aligned/duplicates_removed/{wildcards.sample}.bai {output.bai} &&
            echo 'duplicates_removed bam number of mapped reads:' >> {log} 2>&1 &&
            samtools view -f 2 -c {output.bam} >> {log} 2>&1
            """

else:

    rule handle_duplicates:
        input:
            bam="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam",
            bai="seqnado_output/aligned/blacklist_regions_removed/{sample}.bam.bai",
        output:
            bam=temp("seqnado_output/aligned/duplicates_removed/{sample}.bam"),
            bai=temp("seqnado_output/aligned/duplicates_removed/{sample}.bam.bai"),
        threads: 8
        resources:
            mem="500MB",
        log:
            "seqnado_output/logs/duplicates/{sample}.log",
        script:
            "../scripts/remove_duplicates.py"


if config["shift_atac_reads"]:

    rule shift_atac_alignments:
        input:
            bam="seqnado_output/aligned/duplicates_removed/{sample}.bam",
            bai="seqnado_output/aligned/duplicates_removed/{sample}.bam.bai",
        output:
            bam=temp("seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam"),
            bai=temp(
                "seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam.bai"
            ),
            tmp=temp(
                "seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam.tmp"
            ),
        resources:
            mem="2.5GB",
        threads: 1
        log:
            "seqnado_output/logs/atac_shift/{sample}.log",
        shell:
            """
            rsbamtk shift -b {input.bam} -o {output.tmp} &&
            samtools sort {output.tmp} -@ {threads} -o {output.bam} &&
            samtools index {output.bam} &&
            echo 'Shifted reads' > {log} 2>&1 &&
            samtools view -f 2 -c {output.bam} >> {log} 2>&1
            """

else:

    rule move_bam_to_temp_location:
        input:
            bam="seqnado_output/aligned/duplicates_removed/{sample}.bam",
            bai="seqnado_output/aligned/duplicates_removed/{sample}.bam.bai",
        output:
            bam=temp("seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam"),
            bai=temp(
                "seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam.bai"
            ),
        threads: 1
        log:
            "seqnado_output/logs/atac_shift/{sample}.log",
        shell:
            """
            echo 'Will not shift reads' > {log} &&
            mv {input.bam} {output.bam} &&
            mv {input.bam}.bai {output.bai} &&
            echo 'Number of reads' >> {log} 2>&1 &&
            samtools view -f 2 -c {output.bam} >> {log} 2>&1
            """


rule move_bam_to_final_location:
    input:
        bam="seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam",
        bai="seqnado_output/aligned/shifted_for_tn5_insertion/{sample}.bam.bai",
    output:
        bam="seqnado_output/aligned/{sample,[A-Za-z\\d\\-_]+}.bam",
        bai="seqnado_output/aligned/{sample,[A-Za-z\\d\\-_]+}.bam.bai",
    log:
        "seqnado_output/logs/move_bam/{sample}.log",
    shell:
        """
        mv {input.bam} {output.bam} &&
        mv {input.bai} {output.bai} &&
        echo "BAM moved to final location" > {log} &&
        echo 'Number of reads' > {log} 2>&1 &&
        samtools view -f 2 -c {output.bam} >> {log} 2>&1
        """


localrules:
    move_bam_to_final_location,
