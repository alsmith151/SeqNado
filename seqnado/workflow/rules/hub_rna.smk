import pathlib
import re
import numpy as np

rule generate_hub_for_rnaseq:
    input:
        bigwig=expand(
            "seqnado_output/bigwigs/{method}/{sample}_{strand}.bigWig",
            method=PILEUP_METHODS,
            sample=SAMPLE_NAMES,
            strand=["plus", "minus"],
        ),
        report="seqnado_output/qc/full_qc_report.html",
    output:
        hub=os.path.join(
            config["ucsc_hub_details"]["directory"],
            f"{config['ucsc_hub_details']['name']}.hub.txt",
        ),
    log:
        log=f"seqnado_output/logs/{config['ucsc_hub_details']['name']}.hub.log",
    script:
        "../scripts/create_hub_transcriptomics.py"


localrules:
    generate_hub_for_rnaseq,
