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
import re


##################
# Pipeline setup #
##################

# Read in parameter file
P.get_parameters("config.yml")


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
CREATE_BIGWIGS = P.PARAMS.get('bigwig_create')
CALL_PEAKS = P.PARAMS.get('peaks_call')
CREATE_HUB = P.PARAMS.get('hub_create')

    
# Ensures that all fastq are named correctly
if not os.path.exists('fastq'):
    os.mkdir('fastq')

fastqs = dict()
for fq in glob.glob('*.fastq*'):
    fq_renamed =  (fq.replace('Input', 'input')
                    .replace('INPUT', 'input')
                    .replace('R1.fastq', '1.fastq')
                    .replace('R2.fastq', '2.fastq')
                    )
    
    fastqs[os.path.abspath(fq)] = os.path.join('fastq', fq_renamed)

for src, dest in fastqs.items():
    if not os.path.exists(dest):
        os.symlink(src, dest)


#############
# Read QC   #
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


@follows(mkdir('trimmed'), mkdir("statistics/trimming/data"))
@collate(
    "fastq/*.fastq*",
    regex(r"fastq/(.*)_R?[12].fastq(?:.gz)?"),
    r"trimmed/\1_1_val_1.fq",
)
def fastq_trim(infiles, outfile):

    """Trim adaptor sequences from fastq files using trim_galore"""

    fq1, fq2 = infiles
    fq1_basename, fq2_basename = os.path.basename(fq1), os.path.basename(fq2)

    outdir = 'trimmed'
    trim_options = P.PARAMS.get("trim_options", '')
    cores = P.PARAMS['pipeline_n_cores'] if int(P.PARAMS['pipeline_n_cores']) <= 8 else '8'

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


###############
# Alignment   #
###############


@follows(mkdir("bam"), mkdir('statistics/alignment'), fastq_trim)
@collate(
    "trimmed/*.fq", regex(r"trimmed/(.*)_[12]_val_[12].fq"), r"bam/\1.bam"
)
def fastq_align(infiles, outfile):
    """
    Aligns fq files.

    Uses bowtie2 before conversion to bam file using Samtools view. 
    Bam file is then sorted and the unsorted bam file is replaced.
    
    """

    fq1, fq2 = infiles
    basename = os.path.basename(outfile).replace(".bam", '')
    sorted_bam = outfile.replace(".bam", "_sorted.bam")

    aligner = P.PARAMS.get("aligner_aligner", "bowtie2")
    aligner_options = P.PARAMS.get("aligner_options", "")
    blacklist = P.PARAMS.get("genome_blacklist", "")

    statement = [
        "%(aligner_aligner)s -x %(aligner_index)s -1 %(fq1)s -2 %(fq2)s %(aligner_options)s |",
        "samtools view - -b > %(outfile)s &&",
        "samtools sort -@ %(pipeline_n_cores)s -o %(sorted_bam)s %(outfile)s",
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
        job_queue=P.PARAMS["pipeline_cluster_queue"],
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
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory=P.PARAMS["pipeline_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )


##############
# Mapping QC #
##############

@follows(fastq_align)
@transform(fastq_align, regex(r'.*/(.*).bam'), r'statistics/alignment/\1.txt')
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


#####################
# Remove duplicates #
#####################


@follows(create_bam_index, mkdir("bam_processed"))
@transform(fastq_align, regex(r"bam/(.*.bam)"), r"bam_processed/\1")
def alignments_filter(infile, outfile):
    """Remove duplicate fragments from bam file."""

    alignments_deduplicate = (
        "--ignoreDuplicates" if P.PARAMS.get("alignments_deduplicate") else " "
    )
    alignments_filter_options = P.PARAMS.get("alignments_filter_options")

    breakpoint()
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
            alignments_filter_options if alignments_filter_options else " ",
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
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory=P.PARAMS["pipeline_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )


###########
# BigWigs #
###########

@active_if(CREATE_BIGWIGS)
@follows(mkdir("bigwigs"), alignments_filter)
@transform('bam_processed/*.bam', regex(r"bam_processed/(.*).bam"), r"bigwigs/\1.bigWig")
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
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory=P.PARAMS["pipeline_memory"],
        job_pipeline_n_cores=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


##############
# Call peaks #
##############

@active_if(CALL_PEAKS)
@follows(mkdir("peaks"))
@transform(
    alignments_filter,
    regex(r'.*/(.*?)(?<!input).bam'),
    r"peaks/\1_peaks.narrowPeak",
)
def call_peaks(infile, outfile):

    peaks_options = P.PARAMS.get('peaks_options')
    statement = ['%(peaks_caller)s callpeak -t %(infile)s -n %(outfile)s']

    chipseq_match = re.match(r'.*/(.*)_(.*).bam', infile)

    if chipseq_match:
        samplename = chipseq_match.group(1)
        antibody = chipseq_match.group(2)
        control_file = f'bam_processed/{samplename}_input.bam'
        
        if os.path.exists(control_file):
            statement.append('-c %(input_file)s')
   

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory=P.PARAMS["pipeline_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )


#######################
# UCSC hub generation #
#######################

@transform(call_peaks, regex(r"peaks/(.*).narrowPeak"), r"peaks/\1.bed")
def convert_narrowpeak_to_bed(infile, outfile):

    statement = """awk '{OFS="\\t"; print $1,$2,$3,$4}' %(infile)s > %(outfile)s"""

    P.run(statement, job_queue=P.PARAMS["pipeline_cluster_queue"], job_condaenv=P.PARAMS["conda_env"])


@transform(
    convert_narrowpeak_to_bed,
    regex(r"peaks/(.*).bed"),
    r"peaks/\1.bigBed",
)
def convert_bed_to_bigbed(infile, outfile):

    statement = """bedToBigBed %(infile)s %(genome_chrom_sizes)s %(outfile)s"""
    P.run(statement, job_queue=P.PARAMS["pipeline_cluster_queue"], job_condaenv=P.PARAMS["conda_env"])


@active_if(CREATE_HUB)
@follows(fastq_align, alignments_pileup, alignments_multiqc)
@merge([alignments_pileup, convert_bed_to_bigbed], 
        regex(r'.*'), 
        os.path.join(P.PARAMS.get("hub_dir", ""), P.PARAMS.get("hub_name", "") + ".hub.txt"),
       )
def make_ucsc_hub(infile, outfile, *args):

    import trackhub
    import pickle
    import shutil

    hub_pkl_path =  os.path.join(P.PARAMS['hub_dir'], '.hub.pkl') 

    if os.path.exists(hub_pkl_path) and P.PARAMS.get('hub_append'):

        # Extract previous hub data
        with open(hub_pkl_path, 'rb') as pkl:
            hub, genomes_file, genome, trackdb = pickle.load(pkl)
        
        # Delete previously staged hub

        # Delete symlinks and track db
        for fn in glob.glob(os.path.join(P.PARAMS['hub_dir'], P.PARAMS['genome_name'] + '*')):
            os.unlink(fn)
        
        # Delete hub metadata
        os.unlink(os.path.join(P.PARAMS['hub_dir'], P.PARAMS['hub_name'] + '.hub.txt'))
        os.unlink(os.path.join(P.PARAMS['hub_dir'], P.PARAMS['hub_name'] + '.genomes.txt'))
            
    else:

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
    

    # Move hub to public directory
    trackhub.upload.stage_hub(hub, P.PARAMS['hub_dir'])

    # Save pickle file with data
    with open(hub_pkl_path, 'wb') as pkl:
        pickle.dump([hub, genomes_file, genome, trackdb])



if __name__ == "__main__":
    
    if ("-h" in sys.argv or "--help" in sys.argv):  # If --help then just run the pipeline without setup
        sys.exit(P.main(sys.argv))
    
    elif not 'make' in sys.argv:
        sys.exit(P.main(sys.argv))
    
    elif 'make' in sys.argv:
        sys.exit(P.main(sys.argv))
