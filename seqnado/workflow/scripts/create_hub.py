import os
import pandas as pd
import itertools
import numpy as np
from pathlib import Path
import re
import sys
from loguru import logger
import tracknado

# configre logger to write to snakemake log file if available
if "snakemake" in globals():
    logger.remove()
    logger.add(
        snakemake.log[0],
        format="{time} {level} {message}",
        level="INFO",
    )
    logger.add(sys.stderr, format="{time} {level} {message}", level="ERROR")

def get_rna_samplename(path: str):
    p = Path(path)
    return re.split(r"_[plus|minus]", p.name)[0]

try:
    logger.info("Generating UCSC Genome Browser hub...")
    logger.info(f"Assay: {snakemake.params.assay}, Files: {len(snakemake.input.data)}")
    
    # Set up details
    df = pd.DataFrame(
        snakemake.input.data,
        columns=["fn"],
    )
    
    # Resolve relative paths to absolute paths for consistent parsing
    df["fn"] = df["fn"].apply(lambda x: str(Path(x).resolve()))
    
    # Extract metadata BEFORE TrackFiles processing
    assay_str = str(snakemake.params.assay).split('.')[-1]  # Convert Assay.ATAC -> 'ATAC'
    
    if assay_str in ["ChIP", "CUT&TAG"]:
        df[["samplename", "antibody"]] = df["fn"].str.extract(
            r".*/(.*)_(.*)\.(?:bigBed|bigWig)"
        )
        df["method"] = df["fn"].apply(lambda x: x.split("/")[-3])
        df['norm'] = df['fn'].apply(lambda x: x.split("/")[-2])
    
    elif assay_str == "ATAC":
        df["samplename"] = df["fn"].str.extract(r".*/(.*)\.(?:bigBed|bigWig)")
        df["method"] = df["fn"].apply(lambda x: x.split("/")[-3])
        df["norm"] = df["fn"].apply(lambda x: x.split("/")[-2])
    
    elif assay_str == "RNA":
        df["samplename"] = df["fn"].apply(get_rna_samplename)
        df["method"] = df["fn"].apply(lambda x: x.split("/")[-3])
        df["strand"] = np.where(df["fn"].str.contains("_plus.bigWig"), "plus", "minus")
        df["norm"] = df["fn"].apply(lambda x: x.split("/")[-2])
    
    elif assay_str == 'MCC':
        # Regex pattern to extract method, normalisation, sample, viewpoint
        pattern = re.compile(
        r'seqnado_output/(?:bigwigs|peaks)/'
        r'(?P<method>[^/]+)/'
        r'(?:(?P<norm>[^/]+)/)?'
        r'(?P<samplename>.*?)_(?P<viewpoint>[^/.]+)\.(?:bigWig|bigBed)'
    )
        # Extract the method, normalisation, sample, and viewpoint from the file path
        df_meta = df['fn'].str.extract(pattern)
        df = df.join(df_meta)
    
    logger.debug(f"After metadata extraction - columns: {df.columns.tolist()}")
    
    # Check for missing values in critical columns
    if assay_str == "ATAC":
        for col in ['samplename', 'method', 'norm']:
            null_count = df[col].isnull().sum()
            if null_count > 0:
                logger.error(f"Column '{col}' has {null_count} null values!")
    
    # Use the TrackFiles class to deduplicate files (it will preserve existing columns)
    df = tracknado.TrackFiles(files=df, deduplicate=True).files
    
    # Check that the dataframe is not empty i.e. no files were found
    if df.empty:
        raise ValueError("No bigwigs or bigbeds found in the input directory. Please ensure that create_bigwigs has been set to True in the config file.")
    
    # Create hub design
    design = tracknado.TrackDesign.from_design(
        df,
        color_by=snakemake.params.color_by,
        subgroup_by=snakemake.params.subgroup_by
        if any(snakemake.params.subgroup_by)
        else None,
        supergroup_by=snakemake.params.supergroup_by,
        overlay_by=snakemake.params.overlay_by,
    )
    
    # Generate hub files
    outdir = Path(str(snakemake.output.hub)).parent
    hub = tracknado.HubGenerator(
        track_design=design,
        genome=snakemake.params.genome,
        hub_name=snakemake.params.hub_name,
        description_html=Path(snakemake.input.report),
        hub_email=snakemake.params.hub_email,
        custom_genome=snakemake.params.custom_genome,
        genome_twobit=snakemake.params.genome_twobit,
        genome_organism=snakemake.params.genome_organism,
        genome_default_position=snakemake.params.genome_default_position,
        outdir=outdir,
    )
    
    hub.stage_hub()
    design.to_pickle(outdir / ".track_design.pkl")
    
    logger.info(f"✓ Hub files generated successfully in {outdir}")
    logger.info("="*80)

except Exception as e:
    logger.error("="*80)
    logger.error(f"ERROR: Failed to generate UCSC hub: {e}")
    logger.error(f"Exception type: {type(e).__name__}")
    import traceback
    logger.error(f"Traceback:\n{traceback.format_exc()}")
    logger.error("="*80)
    raise