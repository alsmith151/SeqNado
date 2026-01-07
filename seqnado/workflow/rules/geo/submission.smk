from pathlib import Path
from typing import Any, List


assay_for_protocol = ASSAY.name

rule geo_protocol:
    input:
        [str(p) for p in OUTPUT.files if "/geo_submission/" not in str(p)],
    output:
        OUTPUT_DIR + "/geo_submission/data_processing_protocol.txt",
    params:
        assay=assay_for_protocol,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/geo_protocol.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/geo_protocol.tsv",
    message: "Producing data processing protocol for GEO submission",
    script:
        "../../scripts/produce_data_processing_protocol.py"


rule samples_table:
    input:
        OUTPUT_DIR + "/geo_submission/data_processing_protocol.txt",
    output:
        OUTPUT_DIR + "/geo_submission/samples_table.txt",
    params: 
        output=OUTPUT,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/samples_table.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/samples_table.tsv",
    message: "Generating samples table for GEO submission",
    run:
        from seqnado.outputs.files import GEOFiles
        source_files = [str(p) for p in OUTPUT.files if "/geo_submission/" not in str(p)]

        df = GEOFiles(make_geo_submission_files=True,
                      assay=OUTPUT.assay,
                      design=OUTPUT.design_dataframe,
                      sample_names=OUTPUT.sample_names,
                      config=OUTPUT.config,
                      processed_files=source_files
                      ).metadata

        df.to_csv(output[0], sep="\t", index=False)


rule geo_upload_instructions:
    input:
        OUTPUT_DIR + "/geo_submission/samples_table.txt",
    output:
        instructions=OUTPUT_DIR + "/geo_submission/upload_instructions.txt",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/geo_upload_instructions.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/geo_upload_instructions.tsv",
    message: "Creating GEO upload instructions",
    run:
        import importlib.resources
        import seqnado.data

        source = importlib.resources.files(seqnado.data) / 'geo_upload_instructions.txt'
        with open(source, 'r') as f:
            with open(output.instructions, 'w') as f_out:
                f_out.write(f.read())

def get_files_for_symlink(wc: Any = None) -> List[str]:
    """
    Get all files that need to be symlinked for GEO submission
    """
    from seqnado.outputs.files import GEOFiles
    # Exclude geo_submission files to avoid circular dependencies
    source_files = [str(p) for p in OUTPUT.files if "/geo_submission/" not in str(p)]

    geo_files = GEOFiles(make_geo_submission_files=True,
                         assay=OUTPUT.assay,
                         design=OUTPUT.design_dataframe,
                         sample_names=OUTPUT.sample_names,
                         config=OUTPUT.config,
                         processed_files=source_files)

    fastq_dir = Path(OUTPUT_DIR + "/fastqs")
    fastqs = sorted([str(fastq_dir / fn) for fq_pair in geo_files.raw_files.values() for fn in fq_pair])

    if not geo_files.processed_data_files.empty:
        processed_files = [str(p) for p in geo_files.processed_data_files['path'].tolist()]
    else:
        processed_files = []

    return [*fastqs, *processed_files]

def get_symlinked_files(wc: Any = None) -> List[str]:
    """
    Get all files that have been symlinked for GEO submission
    """
    from seqnado.outputs.files import GEOFiles
    outdir = Path(OUTPUT_DIR + "/geo_submission")

    # Exclude geo_submission files to avoid circular dependencies
    source_files = [str(p) for p in OUTPUT.files if "/geo_submission/" not in str(p)]

    geo_files = GEOFiles(make_geo_submission_files=True,
                         assay=OUTPUT.assay,
                         design=OUTPUT.design_dataframe,
                         sample_names=OUTPUT.sample_names,
                         config=OUTPUT.config,
                         processed_files=source_files)

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
        files=temp(get_symlinked_files()),
    params:
        output=OUTPUT,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/geo_symlink.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/geo_symlink.tsv",
    message: "Creating symlinks for GEO submission files",
    run:
        from pathlib import Path
        from seqnado.outputs.files import GEOFiles

        # Exclude geo_submission files to avoid circular dependencies
        source_files = [str(p) for p in OUTPUT.files if "/geo_submission/" not in str(p)]

        geo_files = GEOFiles(make_geo_submission_files=True,
                             assay=OUTPUT.assay,
                             design=OUTPUT.design_dataframe,
                             sample_names=OUTPUT.sample_names,
                             config=OUTPUT.config,
                             processed_files=source_files)

        fastqs = geo_files.raw_files
        processed_files = geo_files.processed_data_files

        src = Path(OUTPUT_DIR + "/fastqs")
        dest = Path(OUTPUT_DIR + "/geo_submission")

        # Create symlinks for raw files
        for sample_name, fastq in fastqs.items():
            for fq in fastq:
                fq = Path(fq)
                src_file = (src / fq.name).absolute().resolve()
                dest_file = dest / fq.name
                if not dest_file.exists():
                    dest_file.symlink_to(src_file)

        # Create symlinks for processed files
        for _, row in processed_files.iterrows():
            src_file = Path(row['path']).absolute().resolve()
            dest_file = dest / row['output_file_name']
            if not dest_file.exists():
                dest_file.symlink_to(src_file)


rule md5sum:
    input:
        files=get_symlinked_files,
    output:
        OUTPUT_DIR + "/geo_submission/md5sums.txt",
    params:
        geo_dir=OUTPUT_DIR + "/geo_submission",
    log: OUTPUT_DIR + "/logs/geo/md5sum.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/md5sum.tsv",
    message: "Generating MD5 checksums for GEO submission files",
    shell:
        """
        cd {params.geo_dir}

        # Get the basename of the files and store in the infiles variable
        infiles=""
        for f in {input.files}
        do
            infiles="$infiles $(basename $f)"
        done

        md5sum $infiles > md5sums.txt
        """

rule geo_md5_table:
    input:
        md5sums=OUTPUT_DIR + "/geo_submission/md5sums.txt",
    output:
        raw=OUTPUT_DIR + "/geo_submission/raw_data_checksums.txt",
        processed=OUTPUT_DIR + "/geo_submission/processed_data_checksums.txt",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/geo_md5_table.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/geo_md5_table.tsv",
    message: "Creating separate MD5 checksum tables for raw and processed data",
    run:
        import pandas as pd 

        df = pd.read_csv(input.md5sums, sep=r"\s+", header=None, names=["md5sum", "file"])
        df["is_raw"] = df["file"].str.contains(".fastq.gz")

        df_raw = df[df["is_raw"]]
        df_processed = df[~df["is_raw"]]

        for outfile, df_sub in zip([output.raw, output.processed], [df_raw, df_processed]):
            df_sub.rename(columns={"file": "file name", "md5sum": "file checksum"}).to_csv(outfile, index=False, sep="\t")

rule move_to_upload:
    input:
        infiles = get_symlinked_files,
        validated=OUTPUT_DIR + "/geo_submission/.validated",
    output:
        outdir = directory(OUTPUT_DIR + f"/geo_submission/{ASSAY.clean_name}")
    params:
        output=OUTPUT,
    log: OUTPUT_DIR + "/logs/geo/move_to_upload.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/move_to_upload.tsv",
    message: "Moving files to final GEO upload directory",
    shell: """
    mkdir -p {output.outdir}
    for f in {input.infiles}
    do
        cp $f {output.outdir}
    done
    """

rule remove_headers_for_security:
    input:
        infiles = get_symlinked_files
    output:
        validated=OUTPUT_DIR + "/geo_submission/.validated",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/remove_headers.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/remove_headers.tsv",
    message: "Removing headers from TSV files for GEO submission",
    run:
        from pathlib import Path

        for fn in input.infiles:
            path = Path(fn)
            if path.suffix == '.tsv':
                dest = path.with_suffix('.no_headers.tsv')
                
                with open(path, 'r') as f:
                    with open(dest, 'w') as f_out:
                        for line in f:
                            if line.startswith('#'):
                                continue
                            
                            f_out.write(line)
                path.unlink()
                dest.rename(path)
        
        Path(OUTPUT_DIR + "/geo_submission/.validated").touch()



localrules:
    geo_symlink,
    geo_md5_table,
    samples_table,
    geo_upload_instructions,
    move_to_upload,
    remove_headers_for_security
