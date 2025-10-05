from pathlib import Path
import re
import numpy as np

rule generate_hub:
    input:
        data=[
            OUTPUTS.bigwig_files,
            OUTPUTS.bigbed_files,
        ],
        report="seqnado_output/seqnado_report.html",
    output:
        hub=OUTPUTS.ucsc_hub.hub_txt,
    log:
        log=f"seqnado_output/logs/{CONFIG.assay_config.ucsc_hub.name}.hub.log",
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    params:
        assay=ASSAY,
        params=CONFIG.assay_config.ucsc_hub,
        has_consensus_peaks=OUTPUTS.has_consensus_peaks,
    script:
        "../scripts/create_hub.py"


localrules:
    generate_hub