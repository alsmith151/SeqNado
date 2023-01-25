import ngs_pipeline.utils as utils

    
# rule index_filtered_bam:
#     input:
#         bam="aligned_and_filtered/{sample}.bam",
#     output:
#         index="aligned_and_filtered/{sample}.bam.bai",
#     threads:
#         1
#     shell:
#         "samtools index {input.bam} -@ {threads}"

rule feature_counts:
    input:
        bam = expand("aligned_and_filtered/{sample}.bam", sample=SAMPLE_NAMES),
        bai = expand("aligned_and_filtered/{sample}.bam.bai", sample=SAMPLE_NAMES),
        filtering = expand("flags/{sample}.filtering.complete.sentinel", sample=SAMPLE_NAMES),
        annotation = config["genome"]["annotation"]
    output:
        counts = f"feature_counts/read_counts.tsv",
    params:
        options = utils.check_options(config["featurecounts"]["options"]),
    threads:
        config["featurecounts"]["threads"],
    log:
        "logs/readcounts/featurecounts/featurecounts.log"
    shell:
        """
        featureCounts \
        -a \
        {input.annotation} \
        -T \
        {threads} \
        {params.options} \
        -o \
        {output.counts} \
        {input.bam} \
        > {log} 2>&1
        """
