#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:37:44 2019

@author: asmith
"""
import os
import sys
import numpy as np
import pandas as pd
from ruffus import *
from cgatcore import pipeline as P
from ruffus import (
    mkdir,
    follows,
    transform,
    merge,
    originate,
    collate,
    regex,
    add_inputs,
    active_if,
)
from cgatcore.iotools import zap_file, touch_file
from utils import is_none, is_on

##################
# Pipeline setup #
##################

# Read in parameter file
P.get_parameters("config_rna.yml")


# Small edits to config to enable cluster usage
P.PARAMS["cluster_queue_manager"] = P.PARAMS.get("pipeline_cluster_queue_manager")
P.PARAMS["conda_env"] = os.path.basename(os.environ["CONDA_PREFIX"])

# Make sure that params dict is typed correctly
for key in P.PARAMS:
    if is_none(P.PARAMS[key]):
        P.PARAMS[key] = None
    elif is_on(P.PARAMS):
        P.PARAMS[key] = True

# Global variables
CREATE_BIGWIGS = P.PARAMS.get("run_options_bigwigs")
CREATE_HUB = P.PARAMS.get("run_options_hub")

#############
# Pipeline  #
#############


@follows(mkdir("statistics"), mkdir("statistics/fastqc"))
@transform("*.fastq.gz", regex(r"(.*).fastq.gz"), r"statistics/fastqc/\1_fastqc.zip")
def qc_reads(infile, outfile):
    """Quality control of raw sequencing reads"""

    statement = "fastqc -q -t %(pipeline_n_cores)s --nogroup %(infile)s --outdir statistics/fastqc"

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_pipeline_n_cores=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@merge(qc_reads, "statistics/readqc_report.html")
def multiqc_reads(infile, outfile):
    """Collate fastqc reports into single report using multiqc"""

    statement = """export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc statistics/fastqc/ -o statistics -n readqc_report.html"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="2G",
        job_condaenv=P.PARAMS["conda_env"],
    )


######################
# Fastq processing   #
######################


@follows(mkdir("trimmed"), mkdir("statistics/trimming/data"))
@collate(
    "*.fastq.gz",
    regex(r"(.*)_(R)?[12].fastq(?:.gz)?"),
    r"trimming/\1.completed",
)
def fastq_trim(infiles, outfile):
    """Trim adaptor sequences from fastq files using trim_galore"""

    fq1, fq2 = infiles
    fq1_basename, fq2_basename = os.path.basename(fq1), os.path.basename(fq2)

    outdir = "trimmed"
    trim_options = P.PARAMS.get("trim_options", "")
    cores = (
        P.PARAMS["pipeline_n_cores"] if int(P.PARAMS["pipeline_n_cores"]) <= 8 else "8"
    )

    statement = """trim_galore
                   --cores %(cores)s
                   --paired %(trim_options)s
                   --dont_gzip
                   -o %(outdir)s
                   %(fq1)s
                   %(fq2)s
                   """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_pipeline_n_cores=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )
    touch_file(outfile)



###############
# Alignment   #
###############


@follows(mkdir("bam"), mkdir("statistics/alignment"), fastq_trim)
@collate("trimmed/*.fq", regex(r"trimmed/(.*)_R?[12]_val_[12].fq"), r"bam/\1.bam")
def fastq_align(infiles, outfile):
    """
    Aligns fq files.

    Uses STAR before conversion to bam file using Samtools view.
    Bam file is then sorted and the unsorted bam file is replaced.

    """

    basename = os.path.basename(outfile).replace(".bam", "")
    sorted_bam = outfile.replace(".bam", "_sorted.bam")

    blacklist = P.PARAMS.get("genome_blacklist", "")

    statement_align = [
        "STAR",
        "--genomeDir",
        P.PARAMS["aligner_index"],
        "--readFilesIn",
        ",".join(infiles),
        "--readFilesCommand",
        "cat",
        "--outSAMtype",
        "BAM Unsorted",
        "--runThreadN",
        str(P.PARAMS["pipeline_n_cores"]),
        "--outFileNamePrefix",
        outfile.replace(".bam", ""),
        P.PARAMS["aligner_options"] or "",
    ]
    statement_samtools = [
        "samtools",
        "sort",
        "-@",
        str(P.PARAMS["pipeline_n_cores"]),
        "-o",
        sorted_bam,
        f'{outfile.replace(".bam", "")}Aligned.out.bam',
    ]

    if blacklist:
        statement_blacklist = [
            "bedtools",
            "intersect",
            "-v",
            "-a",
            sorted_bam,
            "-b",
            P.PARAMS["blacklist"],
            ">",
            outfile,
            "&&",
            "rm",
            "-f",
            sorted_bam,
        ]
    else:
        statement_blacklist = ["mv", sorted_bam, outfile]

    P.run(
        f'{" ".join(statement_align)} && {" ".join(statement_samtools)} && {" ".join(statement_blacklist)}',
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_pipeline_n_cores=P.PARAMS["pipeline_n_cores"],
        job_memory="32G",
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Zeros the trimmed fastq files
    for fn in infiles:
        zap_file(fn)


@transform(fastq_align, regex(r"bam/(.*)"), r"bam/\1.bai")
def create_bam_index(infile, outfile):
    """Creates an index for the bam file"""

    statement = "samtools index %(infile)s"

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory=P.PARAMS["pipeline_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )


##############
# Mapping QC #
##############


@transform(fastq_align, regex(r".*/(.*).bam"), r"statistics/alignment/\1.txt")
def alignment_statistics(infile, outfile):

    statement = """samtools stats %(infile)s > %(outfile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="2G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(fastq_align, multiqc_reads, alignment_statistics)
@originate("statistics/mapping_report.html")
def alignments_multiqc(outfile):
    """Combines mapping metrics using multiqc"""

    statement = """export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc statistics/alignment/ -o statistics -n alignmentqc_report.html"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="2G",
        job_condaenv=P.PARAMS["conda_env"],
    )


################
# Count reads  #
################


@follows(mkdir("featureCounts"), fastq_align)
@merge(
    fastq_align,
    f'featureCounts/{P.PARAMS["featurecounts_output"]}',
)
def count_reads(infiles, outfile):

    fnames = " ".join(infiles)

    statement = """featureCounts 
                   -a %(featurecounts_gtf)s 
                   -o %(outfile)s 
                   -T %(pipeline_n_cores)s  
                   %(featurecounts_options)s 
                   %(fnames)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_pipeline_n_cores=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


###########
# Bigwigs #
###########


@follows(mkdir("bigwigs"), create_bam_index)
@transform(
    fastq_align,
    regex(r"bam/(.*).bam"),
    r"bigwigs/\1_plus.bigWig",
)
def alignments_pileup(infile, outfile):

    plus = outfile
    minus = outfile.replace("plus", "minus")

    statement = """bamCoverage -b %(infile)s -o %(plus)s --filterRNAstrand forward -p %(pipeline_n_cores)s &&
                   bamCoverage -b %(infile)s -o %(minus)s --filterRNAstrand reverse -p %(pipeline_n_cores)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_pipeline_n_cores=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@active_if(CREATE_HUB)
@follows(
    fastq_align,
    alignments_pileup,
    alignments_multiqc,
)
@merge(
    alignments_pileup,
    regex(r".*"),
    os.path.join(
        P.PARAMS.get("hub_dir", ""), P.PARAMS.get("hub_name", "") + ".hub.txt"
    ),
)
def make_ucsc_hub(infile, outfile, *args):

    import trackhub
    import shutil

    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name=P.PARAMS["hub_name"],
        short_label=P.PARAMS.get("hub_short"),
        long_label=P.PARAMS.get("hub_long"),
        email=P.PARAMS["hub_email"],
        genome=P.PARAMS["genome_name"],
    )

    bigwigs = [fn for fn in infile if ".bigWig" in fn]
    bigbeds = [fn for fn in infile if ".bigBed" in fn]

    for bw in bigwigs:

        track = trackhub.Track(
            name=os.path.basename(bw).replace(".bigWig", ""),
            source=bw,  # filename to build this track from
            visibility="full",  # shows the full signal
            color="128,0,5",  # brick red
            autoScale="on",  # allow the track to autoscale
            tracktype="bigWig",  # required when making a track
        )

        trackdb.add_tracks(track)

    for bb in bigbeds:
        track = trackhub.Track(
            name=os.path.basename(bb).replace(".bigBed", ""),
            source=bb,  # filename to build this track from
            color="0,0,0",  # brick red
            tracktype="bigBed",  # required when making a track
        )

        trackdb.add_tracks(track)

    # Stage the hub
    trackhub.upload.stage_hub(hub=hub, staging="hub_tmp_dir")

    # Copy to the new location
    shutil.copytree(
        "hub_tmp_dir",
        P.PARAMS["hub_dir"],
        dirs_exist_ok=True,
        symlinks=P.PARAMS.get("hub_symlink", False),
    )

    # Delete the staged hub
    shutil.rmtree("hub_tmp_dir")


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
