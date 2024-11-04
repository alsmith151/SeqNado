def get_files_for_symlink(wc):

    fastqs = DESIGN.fastq_paths
    bigwigs = [fn for fn OUTPUT.bigwigs if 'deeptools' in fn]
    peaks = OUTPUT.peaks
    return [*fastqs, *bigwigs, *peaks]


rule geo_symlink:
    input:
        files=get_files_for_symlink,
    output:
        linked_files=direcory("seqnado_output/geo_submission/"),
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
        files=geo_symlink.output.linked_files,
    output:
        "seqnado_output/geo_submission/md5sums.txt",
    shell:
        """
        cd seqnado_output/geo_submission
        md5sum > md5sums.txt
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

        df = pd.read_csv("seqnado_output/geo_submission/md5sums.txt", sep=" ", header=None)
        df.columns = ["md5sum", "file"]

        df = df.assign(is_raw=lambda df: df.file.str.contains(".fastq.gz"))

        df_raw = df[df.is_raw]
        df_processed = df[~df.is_raw]

        for df in ([output.raw, output.processed], [df_raw, df_processed]):
            df = df.rename(columns={"file": "file name", "md5sum": "file checksum"})[['file name', 'file checksum']]
            df.to_csv(output.raw, index=False, sep="\t")


rule samples_table:
    input:
        design=config["design"],
    output:
        "seqnado_output/geo_submission/samples_table.txt",
    params: 
        assay=ASSAY
        pipeline_config=config
        design=DESIGN
    container: None
    run:
        
        df = design.to_geo_dataframe(assay, pipeline_config)
        df.to_csv("seqnado_output/geo_submission/samples_table.txt", sep="\t", index=False)


        

        
        