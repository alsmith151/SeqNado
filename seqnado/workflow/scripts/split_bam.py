import logging
import os
import pysam
import shutil
import subprocess
import sys
from optparse import OptionParser
from loguru import logger

# Set up logging
logger.add(snakemake.log[0], level="INFO")

    __version__ = "1.0.5"

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


    def create_report(output_prefix, stats):
        """
        Create a report file with statistics from the BAM processing in TSV format.

        Parameters:
        output_prefix (str): Prefix used for output files.
        stats (dict): A dictionary containing the statistics to report.
        """
        report_file = output_prefix + "_report.tsv"
        with open(report_file, "w") as report:
            # Writing the headers
            headers = stats.keys()
            report.write("\t".join(headers) + "\n")
            
            # Writing the values
            values = [str(stats[key]) for key in headers]
            report.write("\t".join(values) + "\n")

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


    def main():
        parser = OptionParser(usage="%prog [options]", version="%prog " + __version__)
        parser.add_option("-i", dest="bam_file", help="BAM file of the composite genome")
        parser.add_option("-o", "--output", dest="out_prefix", help="Output prefix")
        parser.add_option("-g", "--sample-prefix", dest="sample_prefix", default="hg38", help="Prefix for exogenous chromosome IDs")
        parser.add_option("-p", "--exo-prefix", dest="chr_prefix", default="dm6", help="Prefix for exogenous chromosome IDs")
        parser.add_option("-q", "--mapq", dest="map_qual", type="int", default=30, help="Mapping quality threshold")


        (options, args) = parser.parse_args()
        if not (options.bam_file and options.out_prefix):
            parser.print_help()
            sys.exit()

        process_bam(options.bam_file, options.out_prefix, options.chr_prefix, options.sample_prefix, options.map_qual)
        stats = process_bam(options.bam_file, options.out_prefix, options.chr_prefix, options.sample_prefix, options.map_qual)
        create_report(options.out_prefix, stats)

with logger.catch():
    if __name__ == "__main__":
        main()
