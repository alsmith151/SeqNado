from pathlib import Path
import re
import numpy as np
from seqnado import Assay
from seqnado.workflow.helpers.visualise import get_hub_input_files

# Only create hub if enabled AND bigwigs/bigbeds will be created
if CONFIG.assay_config.create_ucsc_hub and (OUTPUT.bigwig_files or OUTPUT.bigbed_files or OUTPUT.sentinel_files):
    
    rule generate_hub:
        input:
            data=lambda wc: get_hub_input_files(wc, OUTPUT, ASSAY, rules),
            report=OUTPUT_DIR + "/seqnado_report.html",
        output:
            hub=OUTPUT.ucsc_hub_files,
        params:
            assay=ASSAY,
            has_consensus_peaks=OUTPUT.has_consensus_peaks,
            genome=getattr(CONFIG.assay_config.ucsc_hub, "genome", None),
            hub_name=getattr(CONFIG.assay_config.ucsc_hub, "name", "SeqnadoHub"),
            hub_email=getattr(CONFIG.assay_config.ucsc_hub, "email", None),
            custom_genome=lambda wc: True if getattr(CONFIG.assay_config.ucsc_hub, "two_bit", None) else False,
            genome_twobit=getattr(CONFIG.assay_config.ucsc_hub, "two_bit", None),
            genome_organism=getattr(CONFIG.assay_config.ucsc_hub, "organism", None),
            genome_default_position=getattr(CONFIG.assay_config.ucsc_hub, "default_position", None),
            color_by=getattr(CONFIG.assay_config.ucsc_hub, "color_by", None),
            subgroup_by=getattr(CONFIG.assay_config.ucsc_hub, "subgroup_by", None),
            supergroup_by=getattr(CONFIG.assay_config.ucsc_hub, "supergroup_by", None),
            overlay_by=getattr(CONFIG.assay_config.ucsc_hub, "overlay_by", None),
        containerized: False
        log: OUTPUT_DIR + f"/logs/{getattr(CONFIG.assay_config.ucsc_hub, 'name', 'SeqnadoHub')}.hub.log",
        benchmark: OUTPUT_DIR + f"/.benchmark/visualise/{getattr(CONFIG.assay_config.ucsc_hub, 'name', 'SeqnadoHub')}_hub.tsv",
        message: f"Generating UCSC Genome Browser hub: {getattr(CONFIG.assay_config.ucsc_hub, 'name', 'SeqnadoHub')}"
        script: "../../scripts/create_hub_with_tracknado.py"


    localrules:
        generate_hub