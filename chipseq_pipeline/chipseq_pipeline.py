#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline to process ChIP-seq data from fastq to bigWig generation.

In order to run correctly, input files are required to be in the format: 

SampleName1_(Input|AntibodyUsed)_(R)1|2.fastq.gz

"""

# import packages
import sys
import os
import seaborn as sns
import glob
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
from cgatcore.iotools import zap_file
from utils import is_none, is_on


##################
# Pipeline setup #
##################

# Read in parameter file
P.get_parameters("chipseq_pipeline.yml")


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
CREATE_BIGWIGS = P.PARAMS.get('bigwigs_create')
CALL_PEAKS = P.PARAMS.get('peaks_call')
CREATE_HUB = P.PARAMS.get('hub_create')



def fastq_format():
    """Ensures that all fastq are named correctly
    """

    if os.path.exists('fastq'):
        os.mkdir('fastq')

    fastqs = dict()
    for fq in glob.glob('*.fastq*'):
        fq_renamed =  (fq.replace('Input', 'input')
                        .replace('INPUT', 'input'))
        
        fastqs[os.path.abspath(fq)] = os.path.join('fastq', fq_renamed)

    for src, dest in fastqs.items():
        if not os.path.exists(dest):
            os.symlink(src, dest)


#############
# Read QC   #
#############


@follows(mkdir("statistics"), mkdir("statistics/fastqc"))
@transform("*.fastq.gz", regex(r"(.*).fastq.gz"), r"fastqc/\1_fastqc.zip")
def qc_reads(infile, outfile):

    """Quality control of raw sequencing reads"""

    statement = "fastqc -q -t %(pipeline_n_cores)s --nogroup %(infile)s --outdir statstics/fastqc"

    P.run(
        statement,
        job_queue=P.PARAMS["queue"],
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
        job_queue=P.PARAMS["queue"],
        job_memory="2G",
        job_condaenv=P.PARAMS["conda_env"],
    )


######################
# Fastq processing   #
######################


@follows(mkdir("statistics/trimming/data"))
@collate(
    "fastq/*.fastq*",
    regex(r"(.*)_[12].fastq(?:.gz)?"),
    r"ccanalyser_preprocessing/trimmed/\1_1_val_1.fq.gz",
)
def fastq_trim(infiles, outfile):

    """Trim adaptor sequences from fastq files using trim_galore"""

    fq1, fq2 = infiles
    fq1_basename, fq2_basename = os.path.basename(fq1), os.path.basename(fq2)

    outdir = os.path.dirname(outfile)
    trim_options = (
        P.PARAMS["trim_options"] if not is_none(P.PARAMS["trim_options"]) else ""
    )
    statement = """trim_galore
                   --cores %(pipeline_n_cores)s
                   --paired %(trim_options)s
                   --gzip
                   --dont_gzip
                   -o %(outdir)s
                   %(fq1)s
                   %(fq2)s
                   && mv ccanalyser_preprocessing/trimmed/%(fq1_basename)s_trimming_report.txt statistics/trimming/data
                   && mv ccanalyser_preprocessing/trimmed/%(fq2_basename)s_trimming_report.txt statistics/trimming/data
                   """
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_pipeline_n_cores=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


###############
# Alignment   #
###############


@follows(mkdir("bam"), fastq_trim)
@collate(
    "trimmed/*.fq*", regex(r"trimmed/(.*)_[1|2]_val_[1|2].fq(?:.gz)?"), r"bam/\1.bam"
)
def fastq_align(infiles, outfile):

    """Aligns fq files using bowtie2 before conversion to bam file using
    Samtools view. Bam file is then sorted and the unsorted bam file is replaced"""

    fq1, fq2 = infiles
    basename = os.path.basename(outfile).replace(".bam")
    sorted_bam = outfile.replace(".bam", "_sorted.bam")

    aligner = P.PARAMS.get("aligner_aligner", "bowtie2")
    aligner_options = P.PARAMS.get("aligner_options", "")
    blacklist = P.PARAMS.get("genome_blacklist", "")

    statement = [
        "%(aligner_aligner)s %(aligner_index)s -1 %(fq1)s -2 %(fq2)s %(options)s "
        "2> statistics/alignment/%(basename)s.log |",
        "samtools view -b - > %(outfile)s &&",
        "samtools sort -@ %(pipeline_n_cores)s -m 5G -o %(sorted_bam)s %(outfile)s",
    ]

    if blacklist:
        # Uses bedtools intersect to remove blacklisted regions
        statement.append(
            "&& bedtools intersect -v -b %(blacklist)s -a %(sorted_bam)s > %(outfile)s"
        )
        statement.append("&& rm -f %(sorted_bam)s")

    else:
        statement.append("&& mv %(sorted_bam)s %(outfile)s")

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["queue"],
        job_pipeline_n_cores=P.PARAMS["pipeline_n_cores"],
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
        job_queue=P.PARAMS["queue"],
        job_memory=P.PARAMS["memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )


##############
# Mapping QC #
##############


@originate("report/mapping_report.html")
def alignments_multiqc(infile, outfile):

    """Combines mapping metrics using multiqc"""

    statement = """export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc statistics/alignment/ -o report -n alignment_report.html"""
    P.run(
        statement,
        job_queue=P.PARAMS["queue"],
        job_memory="2G",
        job_condaenv=P.PARAMS["conda_env"],
    )


#####################
# Remove duplicates #
#####################


@follows(create_bam_index, mkdir("bam_processed"))
@transform(fastq_align, regex(r"bam/(.*.bam)"), r"bam_processed/\1")
def alignments_filter(infile, outfile):
    """Remove duplicate fragments from bam file."""

    alignments_deduplicate = (
        "--ignoreDuplicates" if P.PARAMS.get("alignments_deduplicate") else ""
    )
    alignments_filter_options = P.PARAMS.get("alignment_filter_options", "")

    if alignments_deduplicate or alignments_filter_options:

        statement = [
            "alignmentSieve",
            "-b",
            infile,
            "-o",
            outfile,
            "-p",
            "%(pipeline_n_cores)s",
            alignments_deduplicate,
            alignments_filter_options,
            "&& samtools sort -o %(outfile)s.tmp %(outfile)s -@ %(pipeline_n_cores)s",
            "&& mv %(outfile)s.tmp %(outfile)s",
            "&& samtools index %(outfile)s",
            "&& rm -f %(outfile)s.tmp",
        ]

    else:
        infile_abspath = os.path.abspath(infile)
        statement = [
            "ln -s %(infile_abspath)s %(outfile)s",
            "&& ln -s %(infile_abspath)s.bai %(outfile)s.bai",
        ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["queue"],
        job_memory=P.PARAMS["memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )


###########
# BigWigs #
###########


@follows(mkdir("bigwigs"))
@transform(alignments_filter, regex(r"deduplicated/(.*).bam"), r"bigwigs/\1.bigWig")
def alignments_pileup(infile, outfile):

    cmd = [
        "bamCoverage",
        "-b",
        infile,
        "-o",
        outfile,
        "-p",
        "%(pipeline_n_cores)s",
        "%(bigwig_options)s",
    ]

    statement = " ".join(cmd)

    P.run(
        statement,
        job_queue=P.PARAMS["queue"],
        job_memory=P.PARAMS["memory"],
        job_pipeline_n_cores=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


##############
# Call peaks #
##############

@active_if(P.PARAMS.get('peaks_call'))
@follows(mkdir("peaks"))
@transform(
    alignments_filter,
    regex(r"bam_processed/(.*)_(?!input)(.*).bam"),
    r"peaks/\1_\2_peaks.narrowPeak",
    extras=[r'\1', r'\2']
)
def call_peaks(infile, outfile, samplename, antibody):

    peaks_options = P.PARAMS.get('peaks_options')
    input_file = f'bam_processed/{samplename}_input.bam'

    statement = ['%(peaks_caller)s callpeak -t %(infile)s -n %(samplename)s_%(antibody)s --outdir peaks/']

    if os.path.exists(input_file):
        statement.append('-c %(input_file)s')

    P.run(
        statement,
        job_queue=P.PARAMS["queue"],
        job_memory=P.PARAMS["memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )


#######################
# UCSC hub generation #
#######################

@transform(call_peaks, regex(r"peaks/(.*).narrowPeak"), r"peaks/\1.bed")
def convert_narrowpeak_to_bed(infile, outfile):

    statement = """awk '{OFS="\\t"; print $1,$2,$3,$4}' %(infile)s > %(outfile)s"""

    P.run(statement, job_queue=P.PARAMS["queue"], job_condaenv=P.PARAMS["conda_env"])


@transform(
    convert_narrowpeak_to_bed,
    regex(r"peaks/(.*).bed"),
    r"peaks/\1.bigBed",
)
def convert_bed_to_bigbed(infile, outfile):

    statement = """bedToBigBed %(infile)s %(genome_chrom_sizes)s %(outfile)s"""
    P.run(statement, job_queue=P.PARAMS["queue"], job_condaenv=P.PARAMS["conda_env"])



@merge([alignments_pileup, convert_bed_to_bigbed], 
           regex(r'(.*).(?:bigWig|bigBed)'), 
           os.path.join(P.PARAMS.get("hub_dir", ""), P.PARAMS.get("hub_name", "") + ".hub.txt"),
       )
def make_ucsc_hub(infile, outfile):

    import trackhub

    hub, genomes_file, genome, trackdb = trackhub.default_hub(
            hub_name=P.PARAMS["hub_name"],
            short_label=P.PARAMS.get("hub_short"),
            long_label=P.PARAMS.get("hub_long"),
            email=P.PARAMS["hub_email"],
            genome=P.PARAMS["genome_name"],
        )
    
    bigwigs = [fn for fn in infile if '.bigWig' in fn]
    bigbeds = [fn for fn in infile if '.bigBed' in fn]

    for bw in bigwigs:
        
        track = trackhub.Track(
        name=os.path.basename(bw).replace('.bigWig', ''), 
        source=bw,      # filename to build this track from
        visibility='full',  # shows the full signal
        color='128,0,5',    # brick red
        autoScale='on',     # allow the track to autoscale
        tracktype='bigWig', # required when making a track
        )

        trackdb.add_tracks(track)
    
    for bb in bigbeds:
        track = trackhub.Track(
        name=os.path.basename(bb).replace('.bigBed', ''), 
        source=bb,      # filename to build this track from
        color='0,0,0',    # brick red
        tracktype='bigBed', # required when making a track
        )

        trackdb.add_tracks(track)
    

    trackhub.upload.stage_hub(hub, P.PARAMS['hub_dir'])



if __name__ == "__main__":
    
    if ("-h" in sys.argv or "--help" in sys.argv):  # If --help then just run the pipeline without setup
        P.main(sys.argv)
    
    elif not 'make' in sys.argv:
        P.main(sys.argv)
    
    elif 'make' in sys.argv:
        fastq_format()
        P.main(sys.argv)
