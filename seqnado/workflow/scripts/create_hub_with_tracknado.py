import re
import sys
from pathlib import Path
import tracknado as tn
from loguru import logger


def determine_seqnado_assay(parts_lower: list[str]) -> str:
    for i, part in enumerate(parts_lower):
        if part == "seqnado_output":
            return parts_lower[i + 1]

    raise ValueError("Could not determine assay from path")


def from_seqnado_path(path: Path) -> dict[str, str]:
    """
    Extract metadata from seqnado file paths.

    Arguments:
        path (Path): The file path from which to extract metadata.
    Returns:
        dict[str, str]: A dictionary containing extracted metadata fields.

    Assumed that the path follows the seqnado output structure:
        Pattern: .../seqnado_output/{assay}/[bigwigs/peaks]/{method}/{norm}/{sample}_{strand|viewpoint}.[bigWig|bed]
        Example: .../seqnado_output/atac/bigwigs/atac_tn5/cpm/sample1.bigWig

    """
    metadata = {}
    parts = path.parts
    parts_lower = [p.lower() for p in parts]

    metadata["assay"] = determine_seqnado_assay(parts_lower)
    metadata["norm"] = parts[-2]
    metadata["method"] = parts[-3]
    metadata["file_type"] = parts[-4]  # bigwigs or peaks for now

    # samplename is usually the stem, but seqnado sometimes has extensions like .plus/.minus
    # We'll take the first part before any dots or underscores commonly used
    stem = path.stem
    metadata["samplename"] = re.split(r"[._]", stem)[0]

    # If the assay is "MCC" we need to extract the viewpoint from the filename
    # Pattern looks like: /bigwigs/mcc/replicates/{sample}_{viewpoint_group}.bigWig
    if metadata["assay"] == "MCC":
        metadata["viewpoint"] = re.split(r"[._]", stem)[-1].split(".")[0]

    # For RNA we need to extract the strandedness from the filename
    # Pattern looks like: /bigwigs/{method}/{norm}/{sample}_{strand}.bigWig
    elif metadata["assay"] == "RNA":
        metadata["strand"] = re.split(r"[._]", stem)[-1].split(".")[0]

    return metadata


def create_hub_with_tracknado():
    # Configure logger to write to snakemake log file if available
    if "snakemake" in globals():
        logger.remove()
        logger.add(
            snakemake.log[0],
            format="{time} {level} {message}",
            level="INFO",
        )
        logger.add(sys.stderr, format="{time} {level} {message}", level="ERROR")

    try:
        logger.info("🌪️  TrackNado: Generating UCSC Genome Browser hub...")
        logger.info(
            f"Assay: {snakemake.params.assay}, Files: {len(snakemake.input.data)}"
        )

        # 1. Initialize the HubBuilder with the input files
        # We use the fluent API to configure everything in a single chain
        supergroup_by = snakemake.params.supergroup_by
        subgroup_by = snakemake.params.subgroup_by
        overlay_by = snakemake.params.overlay_by
        color_by = snakemake.params.color_by

        builder = (
            tn.HubBuilder()
            .add_tracks(snakemake.input.data)
            # 2. Use the built-in seqnado path extractor
            # This automatically handles directory structures like assay/method/norm/sample.bw
            .with_metadata_extractor(from_seqnado_path)
        )

        if supergroup_by:
            builder = builder.group_by(*supergroup_by, as_supertrack=True)

        if subgroup_by:
            builder = builder.group_by(*subgroup_by)

        if overlay_by:
            builder = builder.overlay_by(*overlay_by)

        if color_by:
            builder = builder.color_by(color_by)
        else:
            builder = builder.color_by(
                "samplename"
            )  # Default coloring by method if none provided

        # If we have a custom genome, set it here
        if (
            hasattr(snakemake.params, "custom_genome")
            and snakemake.params.custom_genome
        ):
            # custom_genome=lambda wc: True if CONFIG.assay_config.ucsc_hub.two_bit else False,
            # genome_twobit=CONFIG.assay_config.ucsc_hub.two_bit,
            # genome_organism=CONFIG.assay_config.ucsc_hub.organism,
            # genome_default_position=CONFIG.assay_config.ucsc_hub.default_position,
            builder = builder.with_custom_genome(
                name=snakemake.params.genome,
                twobit_file=snakemake.params.genome_twobit,
                organism=snakemake.params.genome_organism,
                default_position=snakemake.params.genome_default_position,
            )

        # 4. Build and stage the hub
        # This automatically generates hub.txt, genomes.txt, trackDb.txt
        # and saves a 'tracknado_config.json' sidecar for easy merging later.
        outdir = Path(str(snakemake.output.hub)).parent
        hub = builder.build(
            name=snakemake.params.hub_name,
            genome=snakemake.params.genome,
            outdir=outdir,
            hub_email=snakemake.params.hub_email,
            description_html=Path(snakemake.input.report)
            if hasattr(snakemake.input, "report")
            else None,
        )

        hub.stage_hub()

        logger.info(f"✓ Hub successfully generated in {outdir}")
        logger.info("💡 You can merge multiple hubs later using 'tracknado merge'")

    except Exception as e:
        logger.error("=" * 80)
        logger.error(f"❌ ERROR: Failed to generate UCSC hub: {e}")
        import traceback

        logger.error(f"Traceback:\n{traceback.format_exc()}")
        logger.error("=" * 80)
        raise


if __name__ == "__main__":
    create_hub_with_tracknado()
