"""Helper functions for GEO submission workflows."""

from pathlib import Path
from typing import Any, List


def get_files_for_symlink(OUTPUT, OUTPUT_DIR, wc: Any = None) -> List[str]:
    """
    Get all files that need to be symlinked for GEO submission.

    Args:
        OUTPUT: The OUTPUT object containing files and configuration.
        OUTPUT_DIR: The output directory path.
        wc: Optional wildcards (unused, for compatibility).

    Returns:
        List[str]: List of file paths to symlink.
    """
    from seqnado.outputs.files import GEOFiles

    # Exclude geo_submission files to avoid circular dependencies
    source_files = [str(p) for p in OUTPUT.files if "/geo_submission/" not in str(p)]

    geo_files = GEOFiles(
        make_geo_submission_files=True,
        assay=OUTPUT.assay,
        design=OUTPUT.design_dataframe,
        sample_names=OUTPUT.sample_names,
        config=OUTPUT.config,
        processed_files=source_files,
    )

    fastq_dir = Path(OUTPUT_DIR + "/fastqs")
    fastqs = sorted(
        [
            str(fastq_dir / fn)
            for fq_pair in geo_files.raw_files.values()
            for fn in fq_pair
        ]
    )

    if not geo_files.processed_data_files.empty:
        processed_files = [
            str(p) for p in geo_files.processed_data_files["path"].tolist()
        ]
    else:
        processed_files = []

    return [*fastqs, *processed_files]


def get_symlinked_files(OUTPUT, OUTPUT_DIR, wc: Any = None) -> List[str]:
    """
    Get all files that have been symlinked for GEO submission.

    Args:
        OUTPUT: The OUTPUT object containing files and configuration.
        OUTPUT_DIR: The output directory path.
        wc: Optional wildcards (unused, for compatibility).

    Returns:
        List[str]: List of symlinked file paths.
    """
    from seqnado.outputs.files import GEOFiles

    outdir = Path(OUTPUT_DIR + "/geo_submission")

    # Exclude geo_submission files to avoid circular dependencies
    source_files = [str(p) for p in OUTPUT.files if "/geo_submission/" not in str(p)]

    geo_files = GEOFiles(
        make_geo_submission_files=True,
        assay=OUTPUT.assay,
        design=OUTPUT.design_dataframe,
        sample_names=OUTPUT.sample_names,
        config=OUTPUT.config,
        processed_files=source_files,
    )

    fastqs = [str(outdir / fn) for fqs in geo_files.raw_files.values() for fn in fqs]

    if not geo_files.processed_data_files.empty:
        processed_files = [
            str(outdir / fn)
            for fn in geo_files.processed_data_files["output_file_name"].tolist()
        ]
    else:
        processed_files = []

    return [*fastqs, *processed_files]
