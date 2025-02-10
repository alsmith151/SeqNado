import pathlib
from typing import Any, List

def get_files_for_symlink(wc: Any = None) -> List[str]:
    """
    Get all files that need to be symlinked for GEO submission
    """
    from seqnado.design import GEOFiles
    geo_files = GEOFiles(assay=OUTPUT.assay,
                         design=OUTPUT.design_dataframe,
                         sample_names=OUTPUT.sample_names,
                         config=OUTPUT.config,
                         processed_files=[str(p) for p in OUTPUT.files])

    fastq_dir = pathlib.Path("seqnado_output/fastqs")
    fastqs = sorted([str(fastq_dir / fn) for fq_pair in geo_files.raw_files.values() for fn in fq_pair])
    processed_files = [str(p) for p in geo_files.processed_data_files['path'].tolist()]
    return [*fastqs, *processed_files]

def get_symlinked_files(wc: Any = None) -> List[str]:
    """
    Get all files that have been symlinked for GEO submission
    """
    from seqnado.design import GEOFiles
    outdir = pathlib.Path("seqnado_output/geo_submission")

    geo_files = GEOFiles(assay=OUTPUT.assay,
                         design=OUTPUT.design_dataframe,
                         sample_names=OUTPUT.sample_names,
                         config=OUTPUT.config,
                         processed_files=[str(p) for p in OUTPUT.files])

    fastqs = [str(outdir / fn) for fqs in geo_files.raw_files.values() for fn in fqs]

    if not geo_files.processed_data_files.empty:
        processed_files = [str(outdir / fn) for fn in geo_files.processed_data_files['output_file_name'].tolist()]
    else:
        processed_files = []

    return [*fastqs, *processed_files]


rule geo_symlink:
    input:
        files=get_files_for_symlink,
    output:
        files=get_symlinked_files(),
    params:
        output=OUTPUT,
    container: None
    run:
        import pathlib
        from seqnado.design import GEOFiles

        geo_files = GEOFiles(assay=OUTPUT.assay,
                             design=OUTPUT.design_dataframe,
                             sample_names=OUTPUT.sample_names,
                             config=OUTPUT.config,
                             processed_files=[str(p) for p in OUTPUT.files])

        fastqs = geo_files.raw_files
        processed_files = geo_files.processed_data_files

        src = pathlib.Path("seqnado_output/fastqs")
        dest = pathlib.Path("seqnado_output/geo_submission")

        # Create symlinks for raw files
        for sample_name, fastq in fastqs.items():
            for fq in fastq:
                fq = pathlib.Path(fq)
                src_file = (src / fq.name).absolute().resolve()
                dest_file = dest / fq.name
                if not dest_file.exists():
                    dest_file.symlink_to(src_file)

        # Create symlinks for processed files
        for _, row in processed_files.iterrows():
            src_file = pathlib.Path(row['path']).absolute().resolve()
            dest_file = dest / row['output_file_name']
            if not dest_file.exists():
                dest_file.symlink_to(src_file)


rule md5sum:
    input:
        files=get_symlinked_files,
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

        df = pd.read_csv(input.md5sums, sep=r"\s+", header=None, names=["md5sum", "file"])
        df["is_raw"] = df["file"].str.contains(".fastq.gz")

        df_raw = df[df["is_raw"]]
        df_processed = df[~df["is_raw"]]

        for outfile, df_sub in zip([output.raw, output.processed], [df_raw, df_processed]):
            df_sub.rename(columns={"file": "file name", "md5sum": "file checksum"}).to_csv(outfile, index=False, sep="\t")

rule samples_table:
    output:
        "seqnado_output/geo_submission/samples_table.txt",
    params: 
        output=OUTPUT,
    container: None
    run:
        from seqnado.design import GEOFiles
        df = GEOFiles(assay=OUTPUT.assay,
                      design=OUTPUT.design_dataframe,
                      sample_names=OUTPUT.sample_names,
                      config=OUTPUT.config,
                      processed_files=[str(p) for p in OUTPUT.files]
                      ).metadata
        
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
