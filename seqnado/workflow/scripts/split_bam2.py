import os
import pathlib
import pysam
import shutil
import subprocess
from optparse import OptionParser
from loguru import logger

# Set up logging
logger.add(snakemake.log[0], level="INFO")

def create_headers(bamfile, ex_chr_prefix):
    """Create BAM headers for sample and exogenous genomes."""
    bam_header = bamfile.header
    sample_header, exo_header = {}, {}
    sample_header.update(bam_header)
    exo_header.update(bam_header)

    sample_header["SQ"] = [sq for sq in bam_header["SQ"] if sq["SN"].startswith("chr")]
    exo_header["SQ"] = [sq for sq in bam_header["SQ"] if sq["SN"].startswith(ex_chr_prefix)]

    for header in [sample_header, exo_header]:
        header.setdefault("CO", []).extend([])

    return sample_header, exo_header



def process_bam(bam_file, output_prefix, ex_chr_prefix, sample_genome, map_qual_threshold):
    """Process the BAM file and collect statistics."""
    stats = {
        "bam_file": os.path.basename(bam_file),
        sample_genome + "_reads": 0,
        ex_chr_prefix + "_reads": 0,
        "unmapped_reads": 0,
        "qcfail_reads": 0,
        "duplicate_reads": 0,
        "secondary_reads": 0,
        "low_mapq_reads": 0,
    }

    samfile = pysam.AlignmentFile(bam_file, "rb")
    sample_header, ex_header = create_headers(samfile, ex_chr_prefix)

    with pysam.AlignmentFile(output_prefix + "_" + sample_genome +".bam", "wb", header=sample_header) as sample_out, \
        pysam.AlignmentFile(output_prefix + "_" + ex_chr_prefix + ".bam", "wb", header=ex_header) as exo_out:

        for read in samfile:
            if read.is_unmapped:
                stats["unmapped_reads"] += 1
            elif read.is_qcfail:
                stats["qcfail_reads"] += 1
            elif read.is_duplicate:
                stats["duplicate_reads"] += 1
            elif read.is_secondary:
                stats["secondary_reads"] += 1
            elif read.mapq < map_qual_threshold:
                stats["low_mapq_reads"] += 1
            elif read.reference_name.startswith(ex_chr_prefix):
                stats[ex_chr_prefix + "_reads"] += 1
                exo_out.write(read)
            else:
                stats[sample_genome + "_reads"] += 1
                sample_out.write(read)

    return stats

    report_file = output_prefix + "_report.tsv"
    with open(report_file, "w") as report:
        # Writing the headers
        headers = stats.keys()
        report.write("\t".join(headers) + "\n")
        
        # Writing the values
        values = [str(stats[key]) for key in headers]
        report.write("\t".join(values) + "\n")


with logger.catch():
    logger.info("Split bam files")
    
    bam_file = snakemake.input.bam
    out_prefix = snakemake.params.prefix
    sample_prefix = snakemake.params.genome_prefix
    chr_prefix = snakemake.params.exo_prefix
    map_qual = snakemake.params.map_qual

    process_bam(bam_file, out_prefix, chr_prefix, sample_prefix, map_qual)
