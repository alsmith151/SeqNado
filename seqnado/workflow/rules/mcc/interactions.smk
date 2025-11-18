
rule identify_ligation_junctions:
    input:
        bam="seqnado_output/mcc/{group}/{group}.bam",
        bai="seqnado_output/mcc/{group}/{group}.bam.bai",
    output:
        pairs=temp(expand("seqnado_output/mcc/{{group}}/ligation_junctions/raw/{viewpoint}.pairs", viewpoint=GROUPED_VIEWPOINT_OLIGOS)),
    log:
        "seqnado_output/logs/ligation_junctions/{group}.log",
    threads: 1
    resources:
        mem="1GB",
    container: 'oras://ghcr.io/alsmith151/seqnado_pipeline:latest',
    params:
        outdir="seqnado_output/mcc/{group}/ligation_junctions/raw/",
    shell:
        """
        mccnado identify-ligation-junctions \
        {input.bam} \
        {params.outdir}
        """


rule sort_ligation_junctions:
    input:
        pairs="seqnado_output/mcc/{group}/ligation_junctions/raw/{viewpoint}.pairs",
    output:
        pairs=temp("seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.pairs"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/sort_ligation_junctions/{group}_{viewpoint}.log",
    shell:
        """
        sort -k2,2 -k4,4 -k3,3n -k5,5n {input.pairs} > {output.pairs}
        """

rule bgzip_pairs:
    input:
        pairs="seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.pairs",
    output:
        pairs="seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.pairs.gz",
    log:
        "seqnado_output/logs/bgzip_pairs/{group}_{viewpoint}.log",
    shell:
        """
        bgzip -c {input.pairs} > {output.pairs}
        """

rule make_cooler:
    input:
        pairs="seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.pairs.gz",
        bins="seqnado_output/resources/genomic_bins.bed",
    output:
        cooler=temp("seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.cool"),
    log:
        "seqnado_output/logs/make_cooler/{group}_{viewpoint}.log",
    params:
        resolution=CONFIG.assay_config.mcc.resolution,
        genome=CONFIG.genome.name,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
    shell:
        """
        cooler cload pairs \
        {input.bins} \
        {input.pairs} \
        {output.cooler} \
        --assembly {params.genome} \
        -c1 2 -p1 3 -c2 4 -p2 5 > {log} 2>&1
        """

rule zoomify_cooler:
    input:
        cooler="seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.cool",
    output:
        cooler=temp("seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.mcool"),
    log:
        "seqnado_output/logs/zoomify_cooler/{group}_{viewpoint}.log",
    params:
        resolutions=CONFIG.assay_config.mcc.resolutions,
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        """
        cooler zoomify {input.cooler} -r {params.resolutions} -o {output.cooler} > {log} 2>&1
        """


rule aggregate_coolers:
    input:
        mcools=expand("seqnado_output/mcc/{group}/ligation_junctions/{viewpoint}.mcool", 
                    group=SAMPLE_GROUPINGS.groupings.keys(), 
                    viewpoint=GROUPED_VIEWPOINT_OLIGOS),
    output:
        mcool="seqnado_output/mcc/{group}/{group}.mcool",
    log:
        "seqnado_output/logs/{group}_aggregate_coolers.log",
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    shell:
        """
        mccnado combine-ligation-junction-coolers \
        {input.mcools} \
        {output.mcool}
        """
