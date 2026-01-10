assay_for_protocol = ASSAY.name

protocol_inputs = [
    str(p)
    for p in OUTPUT.files
    if "/geo_submission/" not in str(p) and not str(p).endswith("/protocol.txt")
]

rule protocol:
    input:
        protocol_inputs,
    output:
        OUTPUT_DIR + "/protocol.txt",
    params:
        assay=assay_for_protocol,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/geo_protocol.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/geo_protocol.tsv",
    message: "Producing data processing protocol",
    script:
        "../../scripts/produce_data_processing_protocol.py"