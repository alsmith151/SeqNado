import pathlib
from typing import Any, List

def get_files_for_symlink(wc: Any = None) -> List[str]:

    fastqs = DESIGN.fastq_paths
    bigwigs = OUTPUT.bigwigs
    if assay == "RNA":
        bigwigs = []
    else:
        bigwigs = [f for f in bigwigs if ("deeptools") in str(f)]

    
    peaks = [OUTPUT.peaks[0]] if len(OUTPUT.peaks) > 0 else [] # Just take the first peak file need to change this to handle multiple peak files


    return [*fastqs, *bigwigs, *peaks]


rule geo_symlink:
    input:
        files=get_files_for_symlink,
    output:
        linked_files=[f"seqnado_output/geo_submission/{pathlib.Path(f).name}" for f in get_files_for_symlink()],
    shell:
        """
        mkdir -p seqnado_output/geo_submission

        for file in {input.files}
        do
            ln -s $(realpath $file) seqnado_output/geo_submission/$(basename $file)
        done
        """

rule md5sum:
    input:
        files=rules.geo_symlink.output.linked_files,
    output:
        "seqnado_output/geo_submission/md5sums.txt",
    shell:
        """
        cd seqnado_output/geo_submission
        md5sum * > md5sums.txt
        cd ../..
        """


rule geo_md5_table:
    input:
        md5sums="seqnado_output/geo_submission/md5sums.txt",
    output:
        raw="seqnado_output/geo_submission/raw_data_checksums.txt",
        processed="seqnado_output/geo_submission/processed_data_checksums.txt",
    container: None
    run:
        import pandas as pd 
        import numpy as np

        df = pd.read_csv("seqnado_output/geo_submission/md5sums.txt", sep=r"\s+", header=None)
        df.columns = ["md5sum", "file"]

        df = df.assign(is_raw=lambda df: df.file.str.contains(".fastq.gz"))

        df_raw = df[df.is_raw]
        df_processed = df[~df.is_raw]

        for (outfile, df) in zip([output.raw, output.processed], [df_raw, df_processed]):
            df = df.rename(columns={"file": "file name", "md5sum": "file checksum"})[['file name', 'file checksum']]
            df.to_csv(outfile, index=False, sep="\t")


rule samples_table:
    input:
        design=config["design"],
    output:
        "seqnado_output/geo_submission/samples_table.txt",
    params: 
        assay=ASSAY,
        pipeline_config=config,
        design=DESIGN,
    container: None
    run:
        
        df = params.design.to_geo_dataframe(params.assay, params.pipeline_config)
        df.to_csv(output[0], sep="\t", index=False)
    
rule protocol:
    output:
        "seqnado_output/geo_submission/protocol.txt",
    params:
        assay=ASSAY,
    script:
        "../scripts/produce_data_processing_protocol.py"


        
localrules:
    geo_symlink,
    md5sum,
    geo_md5_table,
    samples_table