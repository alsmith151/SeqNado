#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline to process ChIP-seq data from fastq to bigWig generation.

In order to run correctly, input files are required to be in the format: 

SampleName1_(Input|AntibodyUsed)_(R)1|2.fastq.gz

"""

# import packages
from genericpath import exists
from posixpath import dirname
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
from cgatcore.iotools import zap_file, touch_file
from utils import is_none, is_on
import re
import glob


##################
# Pipeline setup #
##################

# Read in parameter file
P.get_parameters(glob.glob("config_*.yml")[0])


# Small edits to config to enable cluster usage
P.PARAMS["cluster_queue_manager"] = P.PARAMS.get(
    "pipeline_cluster_queue_manager")
P.PARAMS["conda_env"] = os.path.basename(os.environ["CONDA_PREFIX"])

# Make sure that params dict is typed correctly
for key in P.PARAMS:
    if is_none(P.PARAMS[key]):
        P.PARAMS[key] = None
    elif is_on(P.PARAMS):
        P.PARAMS[key] = True

# Global variables
CREATE_BIGWIGS = P.PARAMS.get("run_options_bigwigs")
CALL_PEAKS = P.PARAMS.get("run_options_peaks")
CREATE_HUB = P.PARAMS.get("run_options_hub")
USE_HOMER = P.PARAMS.get("homer_use")
USE_DEEPTOOLS = P.PARAMS.get("deeptools_use")
USE_MACS = P.PARAMS.get("macs_use")
USE_LANCEOTRON = P.PARAMS.get("lanceotron_use")

# Ensures that all fastq are named correctly
if not os.path.exists("fastq"):
    os.mkdir("fastq")

fastqs = dict()
for fq in glob.glob("*.fastq*"):
    fq_renamed = (
        fq.replace("Input", "input")
        .replace("INPUT", "input")
        .replace("R1.fastq", "1.fastq")
        .replace("R2.fastq", "2.fastq")
    )

    fastqs[os.path.abspath(fq)] = os.path.join("fastq", fq_renamed)

for src, dest in fastqs.items():
    if not os.path.exists(dest):
        os.symlink(src, dest)


###################
# Setup functions #
###################


def set_up_chromsizes():
    """
    Ensures that genome chromsizes are present.

    If chromsizes are not provided this function attempts to download them from UCSC.
    The P.PARAMS dictionary is updated with the location of the chromsizes.

    """

    assert P.PARAMS.get("genome_name"), "Genome name has not been provided."

    if P.PARAMS.get("genome_chrom_sizes"):
        pass

    elif os.path.exists("chrom_sizes.txt.tmp"):
        P.PARAMS["genome_chrom_sizes"] = "chrom_sizes.txt.tmp"

    else:
        from pybedtools.helpers import get_chromsizes_from_ucsc

        get_chromsizes_from_ucsc(
            P.PARAMS["genome_name"], "chrom_sizes.txt.tmp")
        P.PARAMS["genome_chrom_sizes"] = "chrom_sizes.txt.tmp"


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

@follows(mkdir('trimmed'))
@transform('fastq/*.fastq*',
           # Regex negates any filenames matching the paired pattern
           regex(r'(?!.*_[12])^fastq/(.*).fastq.gz'),
           r'trimmed/\1_trimmed.fq')
def fastq_trim_single(infile, outfile):

    statement = '''trim_galore 
                   --cores 
                   %(pipeline_n_cores)s 
                   --dont_gzip 
                   %(trim_options)s 
                   -o trimmed 
                   %(infile)s'''

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_pipeline_n_cores=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )

@follows(mkdir("trimmed"), mkdir("statistics/trimming/data"))
@collate(
    "fastq/*.fastq.gz",
    regex(r"fastq/(.*)_R?[12].fastq(?:.gz)?"),
    r"trimmed/\1_1_val_1.fq",
)
def fastq_trim_paired(infiles, outfile):
    """Trim adaptor sequences from fastq files using trim_galore"""

    fq1, fq2 = infiles
    fq1_basename, fq2_basename = os.path.basename(fq1), os.path.basename(fq2)

    outdir = "trimmed"
    trim_options = P.PARAMS.get("trim_options", "")
    cores = (
        P.PARAMS["pipeline_n_cores"] if int(
            P.PARAMS["pipeline_n_cores"]) <= 8 else "8"
    )

    statement = """trim_galore
                   --cores %(cores)s
                   --paired 
                   %(trim_options)s
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

@follows(mkdir("bam"), mkdir("statistics/alignment"), fastq_trim_single, fastq_trim_paired)
@transform(fastq_trim_single, regex(r"trimmed/(.*)_trimmed.fq"), r"bam/\1.bam")
def fastq_align_single(infile, outfile):
    """
    Aligns fq files.

    Uses bowtie2 before conversion to bam file using Samtools view.
    Bam file is then sorted and the unsorted bam file is replaced.

    """

    basename = os.path.basename(outfile).replace(".bam", "")
    sorted_bam = outfile.replace(".bam", "_sorted.bam")

    aligner = P.PARAMS.get("aligner_aligner", "bowtie2")
    aligner_options = P.PARAMS.get("aligner_options", "")
    blacklist = P.PARAMS.get("genome_blacklist", "")

    statement = [
        "%(aligner_aligner)s -x %(aligner_index)s -U %(infile)s %(aligner_options)s |",
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
    zap_file(infile)


@follows(mkdir("bam"), mkdir("statistics/alignment"), fastq_trim_paired, fastq_trim_single)
@collate("trimmed/*.fq", regex(r"trimmed/(.*)_[12]_val_[12].fq"), r"bam/\1.bam")
def fastq_align_paired(infiles, outfile):
    """
    Aligns fq files.

    Uses bowtie2 before conversion to bam file using Samtools view.
    Bam file is then sorted and the unsorted bam file is replaced.

    """

    fq1, fq2 = infiles
    basename = os.path.basename(outfile).replace(".bam", "")
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


@transform([fastq_align_single, fastq_align_paired], regex(r"bam/(.*)"), r"bam/\1.bai")
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


@transform([fastq_align_single, fastq_align_paired], regex(r".*/(.*).bam"), r"statistics/alignment/\1.txt")
def alignment_statistics(infile, outfile):

    statement = """samtools stats %(infile)s > %(outfile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="2G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(multiqc_reads, alignment_statistics)
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
@transform([fastq_align_single, fastq_align_paired], regex(r"bam/(.*.bam)"), r"bam_processed/\1")
def alignments_filter(infile, outfile):
    """Remove duplicate fragments from bam file."""

    alignments_deduplicate = P.PARAMS.get("alignments_deduplicate")
    alignments_filter_options = P.PARAMS.get("alignments_filter_options")

    if alignments_deduplicate or alignments_filter_options:

        statement = [
            "alignmentSieve",
            "-b",
            infile,
            "-o",
            outfile,
            "-p",
            "%(pipeline_n_cores)s",
            "--ignoreDuplicates" if alignments_deduplicate else "",
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


@active_if(USE_HOMER)
@follows(mkdir("tag/"))
@transform(alignments_filter, regex(r".*/(.*).bam"), r"tag/\1")
def create_tag_directory(infile, outfile):

    statement = [
        "makeTagDirectory",
        outfile,
        P.PARAMS["homer_tagdir_options"] or " ",
        infile,
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


@active_if(CREATE_BIGWIGS and USE_DEEPTOOLS)
@follows(mkdir("bigwigs/deeptools/"))
@transform(
    alignments_filter, regex(
        r"bam_processed/(.*).bam"), r"bigwigs/deeptools/\1_deeptools.bigWig")
def alignments_pileup_deeptools(infile, outfile):

    statement = [
        "bamCoverage",
        "-b",
        infile,
        "-o",
        outfile,
        "-p",
        "%(pipeline_n_cores)s",
        P.PARAMS.get("deeptools_bamcoverage_options") or " ",
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory=P.PARAMS["pipeline_memory"],
        job_pipeline_n_cores=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("bigwigs/homer/"))
@active_if(CREATE_BIGWIGS and USE_HOMER)
@transform(
    create_tag_directory, regex(r".*/(.*)"), r"bigwigs/homer/\1_homer.bigWig", extras=[r"\1"]
)
def alignments_pileup_homer(infile, outfile, tagdir_name):

    outdir = os.path.dirname(outfile)

    statement_bw = [
        "makeBigWig.pl",
        infile,
        P.PARAMS["genome_name"],
        "-chromSizes",
        P.PARAMS["genome_chrom_sizes"],
        "-url",
        P.PARAMS.get("homer_makebigwig_options") or "INSERT_URL_HERE",
        "-webdir",
        outdir,
    ]

    statement_mv = ["&&", 'mv', outfile.replace('_homer.bigWig', '.ucsc.bigWig'), outfile]

    P.run(
        " ".join([*statement_bw, *statement_mv]),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_pipeline_n_cores=1,
        job_condaenv=P.PARAMS["conda_env"],
    )

    try:
        # Rename bigwigs to remove ucsc
        bigwig_src = os.path.join(outdir, f"{tagdir_name}.ucsc.bigWig")
        bigwig_dest = os.path.join(outdir, f"{tagdir_name}.bigWig")
        os.rename(bigwig_src, bigwig_dest)
    except OSError:
        pass


##############
# Call peaks #
##############


@active_if(CALL_PEAKS and USE_MACS)
@follows(mkdir("peaks/macs"))
@transform(
    alignments_filter,
    regex(r".*/(.*?)(?<!input).bam"),
    r"peaks/macs/\1_peaks.narrowPeak",
)
def call_peaks_macs(infile, outfile):

    output_prefix = outfile.replace("_peaks.narrowPeak", "")
    statement = [
        "%(macs_caller)s",
        "callpeak",
        "-t",
        "%(infile)s",
        "-n",
        "%(output_prefix)s",
        P.PARAMS.get("macs_callpeak_options") or " ",
    ]

    chipseq_match = re.match(r".*/(.*)_(.*).bam", infile)

    if chipseq_match:
        samplename = chipseq_match.group(1)
        antibody = chipseq_match.group(2)
        control_file = f"bam_processed/{samplename}_input.bam"

        if os.path.exists(control_file):
            statement.append(f"-c {control_file}")

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory=P.PARAMS["pipeline_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )

@active_if(CALL_PEAKS and USE_LANCEOTRON)
@follows(mkdir("peaks/macs"))
@transform(
    alignments_pileup_deeptools,
    regex(r".*/(.*?)(?<!input).bigWig"),
    r"peaks/lanceotron/\1_L-tron.bed",
)
def call_peaks_lanceotron(infile, outfile):

    chipseq_match = re.match(r".*/(.*)_(.*).bigWig", infile)

    if chipseq_match:
        samplename = chipseq_match.group(1)
        antibody = chipseq_match.group(2)
        control_file = f"bigwigs/deeptools/{samplename}_input.bigWig"

    statement = [
        "lanceotron",
        "callPeaks" if not os.path.exists(control_file) else "callPeaksInput",
        "%(infile)s",
        "-i %(control_file)s" if os.path.exists(control_file) else " ",
        "-f",
        os.path.dirname(outfile),
        "--format bed",
        P.PARAMS.get("lanceotron_callpeak_options") or " ",
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory=P.PARAMS["pipeline_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@active_if(CALL_PEAKS and USE_HOMER)
@follows(mkdir("peaks/homer"))
@transform(
    create_tag_directory,
    regex(r".*/(?!.*_input)(.*)"),
    r"peaks/homer/\1_homer_peaks.bed",
)
def call_peaks_homer(infile, outfile):

    tmp = outfile.replace(".bed", ".txt")
    statement = [
        "findPeaks",
        infile,
        P.PARAMS["homer_findpeaks_options"] or " ",
        "-o",
        tmp,
    ]

    # Finds the matching input file if one extists
    chipseq_match = re.match(r".*/(.*)_(.*)", infile)
    if chipseq_match:
        samplename = chipseq_match.group(1)
        antibody = chipseq_match.group(2)
        control = f"tag/{samplename}_input"

        if os.path.exists(control):
            statement.append(f"-i {control}")

    # Need to convert homer peak format to bed
    statement.append(f"&& pos2bed.pl {tmp} -o {outfile}")

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory=P.PARAMS["pipeline_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )


#######################
# UCSC hub generation #
#######################


@transform(call_peaks_macs, regex(r"peaks/(.*).narrowPeak"), r"peaks/\1.bed")
def convert_narrowpeak_to_bed(infile, outfile):

    statement = """awk '{OFS="\\t"; print $1,$2,$3,$4}' %(infile)s > %(outfile)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@transform(
    [convert_narrowpeak_to_bed, call_peaks_homer, call_peaks_lanceotron],
    regex(r"(.*)/(.*).bed"),
    r"\1/\2.bigBed",
)
def convert_bed_to_bigbed(infile, outfile):

    statement = """cat %(infile)s 
                   | sort -k1,1 -k2,2n > %(infile)s.tmp 
                   && bedToBigBed %(infile)s.tmp %(genome_chrom_sizes)s %(outfile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@active_if(CREATE_HUB)
@follows(
    fastq_align_single,
    fastq_align_paired,
    alignments_pileup_deeptools,
    alignments_pileup_homer,
    alignments_multiqc,
)
@merge(
    [alignments_pileup_deeptools, alignments_pileup_homer, convert_bed_to_bigbed],
    regex(r".*"),
    os.path.join(
        P.PARAMS.get("hub_dir", ""), P.PARAMS.get("hub_name", "") + ".hub.txt"
    ),
)
def make_ucsc_hub(infile, outfile, *args):

    import trackhub
    import shutil
    import seaborn as sns

    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name=P.PARAMS["hub_name"],
        short_label=P.PARAMS.get("hub_short"),
        long_label=P.PARAMS.get("hub_long"),
        email=P.PARAMS["hub_email"],
        genome=P.PARAMS["genome_name"],
    )

    bigwigs = [fn for fn in infile if ".bigWig" in fn]
    colours = dict(zip(bigwigs, sns.color_palette("hls", len(set(bigwigs)))))
    bigbeds = [fn for fn in infile if ".bigBed" in fn]

    for bw in bigwigs:

        track = trackhub.Track(
            name=os.path.basename(bw).replace(".bigWig", ""),
            source=bw,  # filename to build this track from
            visibility="full",  # shows the full signal
            color=",".join(
                    [str(int(x * 255)) for x in colours[bw]]),  # brick red
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


@follows(
    alignments_pileup_homer,
    alignments_pileup_homer,
    call_peaks_homer,
    call_peaks_macs,
    call_peaks_lanceotron,
)
@originate("pipeline_complete.txt")
def full(outfile):
    touch_file(outfile)


if __name__ == "__main__":

    if (
        "-h" in sys.argv or "--help" in sys.argv
    ):  # If --help then just run the pipeline without setup
        sys.exit(P.main(sys.argv))

    elif not "make" in sys.argv:
        sys.exit(P.main(sys.argv))

    elif "make" in sys.argv:
        set_up_chromsizes()
        sys.exit(P.main(sys.argv))
