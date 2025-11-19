from pathlib import Path
import re
import numpy as np

if CONFIG.assay_config.create_ucsc_hub:

    rule generate_hub:
        input:
            data=[
                OUTPUT.bigwig_files,
                OUTPUT.bigbed_files,
            ],
            report=OUTPUT_DIR + "/seqnado_report.html",
        output:
            hub=OUTPUT.ucsc_hub_files,
        log:
            log=OUTPUT_DIR + "/logs/{CONFIG.assay_config.ucsc_hub.name}.hub.log",
        container:
            "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        params:
            assay=ASSAY,
            params=CONFIG.assay_config.ucsc_hub,
            has_consensus_peaks=OUTPUT.has_consensus_peaks,
        script:
            "../scripts/create_hub.py"


    localrules:
        generate_hub