rule bed_to_bigbed:
    input: 
        bed = "peaks/{directory}/{sample}.bed",
    output:
        bigbed = "peaks/{directory}/{sample}.bigBed"
    params:
        chrom_sizes = config["genome"]["chromosome_sizes"]
    shell:
        """
        sort -k1,1 -k2,2n {input.bed} > {input.bed}.tmp &&
        bedToBigBed {input.bed}.tmp {params.chrom_sizes} {output.bigbed}
        """

rule generate_hub:
    input:
        bigbed = expand("peaks/{method}/{sample}.bigBed", method=PEAK_CALL_METHODS, sample=SAMPLE_NAMES_NO_READ),
        bigwig = expand("bigwigs/{method}/{sample}.bigWig", method=PILEUP_METHODS, sample=SAMPLE_NAMES_NO_READ),
        report = "qc/full_qc_report.html",
    output:
        hub = f"{config['ucsc_hub_details']['directory']}/{config['ucsc_hub_details']['name']}.hub.txt",
    log:
        log = f"logs/{config['ucsc_hub_details']['name']}.hub.log"
    # conda:
    #     "../../environment_ngs.yml"
    run:
        
        import pandas as pd
        import itertools

        df = pd.DataFrame(itertools.chain.from_iterable([files for files in [input.bigbed, input.bigwig]]), columns=["filename"])
        df[["samplename", "antibody"]] = df["filename"].str.extract(r".*/(.*)_(.*)\.(?:bigBed|bigWig)")
        df["method"] = df["filename"].apply(lambda x: x.split("/")[-2])

        file_details = f"{os.path.dirname(output.hub)}/hub_details.tsv"        
        df.set_index("filename").to_csv(file_details, sep="\t")

        cmd = ["make-ucsc-hub", 
               " ".join(df["filename"]),
               "-d",
               file_details, 
               "-o", 
               os.path.dirname(output.hub),
               "--hub-name",
               config['ucsc_hub_details']["name"],
               "--hub-email",
               config['ucsc_hub_details']["email"],
               "--genome-name",
               config["genome"]["name"],
               "--description-html",
                input.report,
        ]

        cmd_log = f"echo {' '.join(cmd)} > {log.log}"

        shell(" ".join(cmd))
        shell(cmd_log)


        



