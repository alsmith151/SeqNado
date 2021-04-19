#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 10:47:03 2020

@author: asmith
"""

# import packages
import sys
import os
import re
import gzip
import seaborn as sns
import glob
from cgatcore import pipeline as P
from ruffus import mkdir, follows, transform, merge, originate, collate, split, regex, add_inputs, suffix, active_if
from cgatcore.iotools import zap_file
import click

# Read in parameter file
P.get_parameters('chipseq_pipeline.yml')

# Global variables
use_lanceotron = True if P.PARAMS['use_lanceotron'] else False
use_macs3 = True if not use_lanceotron else False
hub_dir = os.path.join(P.PARAMS["hub_publoc"], P.PARAMS['hub_name'])
assembly_dir = os.path.join(hub_dir, P.PARAMS['hub_genome'])

# Hack to enable slurm cluster
P.PARAMS['cluster_queue_manager'] = 'slurm'
P.PARAMS["conda_env"] = P.PARAMS.get(
    "conda_env", os.path.basename(os.environ["CONDA_PREFIX"])
)


'''Requires input files to be in the format: .*_[input|.*]_[1|2]'''

@follows(mkdir('fastqc'))
@transform('*.fastq.gz', 
           regex(r'(.*).fastq.gz'), 
           r'fastqc/\1_fastqc.zip')
def qc_reads(infile, outfile):
    '''Quality control of raw sequencing reads'''
    statement = 'fastqc -q -t %(threads)s --nogroup %(infile)s --outdir fastqc'
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'],
          job_condaenv=P.PARAMS["conda_env"])


@follows(mkdir('report'))
@merge(qc_reads, 'report/readqc_report.html')
def multiqc_reads (infile, outfile):
    '''Collate fastqc reports into single report using multiqc'''
    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc fastqc/ -o report -n readqc_report.html'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_memory='16G',
          job_condaenv=P.PARAMS["conda_env"])
    
    
@follows(mkdir('trimmed'))
@transform('*.fastq.gz', 
           regex(r'(?!.*_.*_[12])^(.*).fastq.gz'), # Regex negates any filenames matching the paired pattern
           r'trimmed/\1_trimmed.fq.gz')
def trim_reads_single(infile, outfile):
    
    statement = '''trim_galore --cores %(threads)s %(trim_options)s 
                  -o trimmed %(infile)s'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'],
          job_condaenv=P.PARAMS["conda_env"])
    

@collate(r'*.fastq.gz',
         regex(r'(.*)_[1|2].fastq.gz'), 
         r'trimmed/\1_1_val_1.fq.gz')
def trim_reads_paired(infiles, outfile):
    '''Trim adaptor sequences using Trim-galore in paired end mode'''
    fastq1, fastq2 = infiles
    statement = '''trim_galore --cores %(threads)s --paired %(trim_options)s 
                  -o trimmed %(fastq1)s %(fastq2)s'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'],
          job_condaenv=P.PARAMS["conda_env"])
    

@follows(mkdir('bam'), trim_reads_single)
@transform(trim_reads_single, 
           regex(r'trimmed/(.*)_trimmed.fq.gz'), 
           r'bam/\1.bam')
def align_reads_single(infile, outfile):
    ''' Aligns digested fq files using bowtie2'''
    
    sorted_bam = outfile.replace('.bam', '_sorted.bam')
    blacklist = ''
        
    if 'bowtie2_options' in P.PARAMS.keys():
        options = P.PARAMS['bowtie2_options'] if P.PARAMS['bowtie2_options'] else ' '
    
    if 'genome_blacklist' in P.PARAMS.keys():
        blacklist = P.PARAMS['genome_blacklist']
    
    
    statement = ['bowtie2 -x %(bowtie2_index)s -U %(infile)s -p %(threads)s %(options)s |',
                 'samtools view -b - > %(outfile)s &&',
                 'samtools sort -@ %(threads)s -m 5G -o %(sorted_bam)s %(outfile)s']
    
    if blacklist:
        # Uses bedtools intersect to remove blacklisted regions
        statement.append(''' && bedtools intersect -v -b %(blacklist)s -a %(sorted_bam)s > %(outfile)s &&
                          rm -f %(sorted_bam)s''')
    else:
        statement.append('&& mv %(sorted_bam)s %(outfile)s')
    
    P.run(' '.join(statement), 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'],
          job_condaenv=P.PARAMS["conda_env"])
    
    zap_file(infile)


   
@follows(mkdir('bam'), trim_reads_paired)
@collate('trimmed/*.fq.gz', 
         regex(r'trimmed/(.*)_[1|2]_val_[1|2].fq.gz'), 
         r'bam/\1.bam')
def align_reads_paired(infiles, outfile):
    
    ''' Aligns fq files using bowtie2 before conversion to bam file using
        Samtools view. Bam file is then sorted and the unsorted bam file is replaced'''
    
    fq1, fq2 = infiles   
    sorted_bam = outfile.replace('.bam', '_sorted.bam')
    options = ''
    blacklist= None
        
    if P.PARAMS.get('bowtie2_options') not in [None, 'None', '']:
        options = P.PARAMS['bowtie2_options']
    
    if 'genome_blacklist' in P.PARAMS.keys():
        blacklist = P.PARAMS['genome_blacklist']
    
    
    statement = ['bowtie2 -x %(bowtie2_index)s -1 %(fq1)s -2 %(fq2)s -p %(threads)s %(options)s |',
                 'samtools view -b - > %(outfile)s &&',
                 'samtools sort -@ %(threads)s -m 5G -o %(sorted_bam)s %(outfile)s']
    
    if blacklist:
        # Uses bedtools intersect to remove blacklisted regions
        statement.append('''&& bedtools intersect -v -b %(blacklist)s -a %(sorted_bam)s > %(outfile)s &&
                          rm -f %(sorted_bam)s''')
    else:
        statement.append('&& mv %(sorted_bam)s %(outfile)s')
    
    P.run(' '.join(statement), 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'],
          job_condaenv=P.PARAMS["conda_env"])
    
    for fn in infiles:
        zap_file(fn)


@transform([align_reads_single, align_reads_paired], 
           regex(r'bam/(.*)'),
           r'bam/\1.bai')
def create_bam_file_index(infile, outfile):
    """Bam files are compressed. The index allows fast access to different
    slices of the file."""
    statement = 'samtools index %(infile)s'
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'],
          job_condaenv=P.PARAMS["conda_env"])


@transform([align_reads_single, align_reads_paired], 
           regex(r'bam/(.*).bam'), 
           r'bam/\1.picard.metrics')
def mapping_qc(infile, outfile):
    '''Uses picard CollectAlignmentSummaryMetrics to get mapping information.'''
    
    cmd = ['picard CollectAlignmentSummaryMetrics',
           'R=%(genome_fasta)s I=%(infile)s O=%(outfile)s',
           '&> %(outfile)s.log',
           ]
    
    statement = ' '.join(cmd)
         
    P.run(statement, 
          job_queue=P.PARAMS['queue'],
          job_condaenv=P.PARAMS["conda_env"])

@merge(mapping_qc, 'report/mapping_report.html')
def mapping_multiqc (infile, outfile):
    '''Combines mapping metrics using multiqc'''
    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc bam/ -o report -n mapping_report.html'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_memory='16G',
          job_condaenv=P.PARAMS["conda_env"])


@follows(create_bam_file_index, mkdir('deduplicated'), mkdir('tmp'))    
@transform(align_reads_paired,
           regex(r'bam/(.*)'),
           r'deduplicated/\1')
def remove_duplicates_paired(infile, outfile):
    stats = outfile.replace('.bam', '_marked_duplicates.txt')
    statement = ' '.join(['picard MarkDuplicates',
                          'I=%(infile)s',
                          'O=%(outfile)s',
                          'M=%(stats)s',
                          'REMOVE_DUPLICATES=true',
                          'REMOVE_SEQUENCING_DUPLICATES=true',
                          'CREATE_INDEX=true',
                          'TMP_DIR=tmp/'])
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'],
          job_condaenv=P.PARAMS["conda_env"])

@follows(create_bam_file_index, mkdir('deduplicated'), mkdir('tmp'))    
@transform(align_reads_single,
           regex(r'bam/(.*)'),
           r'deduplicated/\1')
def remove_duplicates_single(infile, outfile):
    stats = outfile.replace('.bam', '_marked_duplicates.txt')
    statement = ' '.join(['picard MarkDuplicates',
                          'I=%(infile)s',
                          'O=%(outfile)s',
                          'M=%(stats)s',
                          'REMOVE_DUPLICATES=true',
                          'REMOVE_SEQUENCING_DUPLICATES=true',
                          'CREATE_INDEX=true',
                          'TMP_DIR=tmp/'])
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'],
          job_condaenv=P.PARAMS["conda_env"])



@follows(mkdir('bigwigs'))
@transform(remove_duplicates_paired, regex(r'deduplicated/(.*).bam'), r'bigwigs/\1.bigWig')
def make_bigwig_paired(infile, outfile):
    
    cmd = [ 'bamCoverage', 
            '-b', infile, 
            '-o', outfile,
            '-p', '%(threads)s',
            '%(bigwig_options)s']
    
    statement = ' '.join(cmd)
    
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'],
          job_threads=P.PARAMS['threads'],
          job_condaenv=P.PARAMS["conda_env"])
    
@follows(mkdir('bigwigs'))
@transform(remove_duplicates_single, regex(r'deduplicated/(.*).bam'), r'bigwigs/\1.bigWig')
def make_bigwig_single(infile, outfile):
    
    cmd = [ 'bamCoverage', 
            '-b', infile, 
            '-o', outfile,
            '-p', '%(threads)s',
            '%(bigwig_options)s']
    
    statement = ' '.join(cmd)
    
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'],
          job_threads=P.PARAMS['threads'],
          job_condaenv=P.PARAMS["conda_env"])   
    


@follows(mkdir('peaks'))
@active_if(use_macs3)
@transform([remove_duplicates_single, remove_duplicates_paired], 
           regex(r'deduplicated/(.*)_(?!input|Input|INPUT)(.*).bam'),
           r'peaks/\1_\2_peaks.narrowPeak')
def call_peaks(infile, outfile):
    
    treatment = infile
    treatment_name = os.path.basename(treatment).split('_')[0]
    file_base = outfile.replace('_peaks.narrowPeak', '')
    macs_options = ' '
    
    if P.PARAMS['macs_options']:
        macs_options = P.PARAMS['macs2_options']

    statement = '''macs3 callpeak %(macs_options)s -g %(genome_size)s 
                  -t %(treatment)s -n %(file_base)s'''
    
    
    control_files = glob.glob(f'deduplicated/{treatment_name}_input.bam')
    if control_files:
        control = control_files[0]
        statement += ' -c %(control)s'
      

    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'],
          job_condaenv=P.PARAMS["conda_env"])



@active_if(use_lanceotron)
@transform([make_bigwig_paired, make_bigwig_single],
           regex(r'bigwigs/(.*)_(?!input|Input|INPUT)(.*).bigWig'),
           r'peaks/\1_\2_lanceotron_peaks.bed')
def call_peaks_lanceotron(infile, outfile):
    
    lanceotron_options = ' '
    if P.PARAMS['lanceotron_options']:
        lanceotron_options = P.PARAMS['lanceotron_options']
    
    
    base_name = os.path.basename(infile).replace('.bigWig', '')
    dir_name = f'peaks/{base_name.replace("_lanceotron_peaks.bed", "")}/'


    # Switching to using a filtered peakset (currently 18.11.20 using all called)
    # Columns of called peaks are:
    header = ['chrom',
             'start',
             'end',
             'H3K4me1_score',
             'noise_score',
             'ATAC_score',
             'H3K4me3_score',
             'TF_score',
             'H3K27ac_score']
    
    tsv_header = ' '.join(header)


    
    statement = '''rm -rf %(dir_name)s/merge/ &&
                   python /t1-data/user/lhentges/lanceotron/lanceotron_genome.py
                   %(lanceotron_options)s -f %(dir_name)s %(infile)s &&
                   cat %(dir_name)s/merge/*.bed
                   '''
    
    if P.PARAMS.get('lanceotron_filter') not in [None, '', 'None']:

        statement += '''| python /home/nuffmed/asmith/Data/Projects/chipseq_pipeline/filter_tsv.py
                        -
                        --col_names %(tsv_header)s
                        -q "%(lanceotron_filter)s"
                        '''

    statement += '| sort -k1,1 -k2,2n | cut -f 1-3 > %(outfile)s'

                   
    
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'],
          job_condaenv=P.PARAMS["conda_env"]
         )


@active_if(use_macs3)
@transform(call_peaks, regex(r'peaks/(.*).narrowPeak'), r'peaks/\1.bed')
def convert_narrowpeak_to_bed(infile, outfile):
    
    statement = '''awk '{OFS="\\t"; print $1,$2,$3,$4}' %(infile)s > %(outfile)s'''
    
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_condaenv=P.PARAMS["conda_env"])
    
@transform([convert_narrowpeak_to_bed, call_peaks_lanceotron],
           regex(r'peaks/(.*).bed'),
           r'peaks/\1.bigBed')    
def convert_bed_to_bigbed(infile, outfile):
    
    statement = '''bedToBigBed %(infile)s %(genome_chrom_sizes)s %(outfile)s'''
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_condaenv=P.PARAMS["conda_env"])
        
@follows(mkdir(hub_dir), 
         mkdir(assembly_dir), 
         make_bigwig_single, 
         make_bigwig_paired, 
         convert_bed_to_bigbed)
@originate(os.path.join(hub_dir, 'hub.txt'))
def generate_hub_metadata(outfile):

    content = {'hub': P.PARAMS['hub_name'],
               'shortLabel': P.PARAMS['hub_short'] if P.PARAMS['hub_short'] else P.PARAMS['hub_name'],
               'longLabel': P.PARAMS['hub_long'] if P.PARAMS['hub_long'] else P.PARAMS['hub_name'],
               'genomesFile': 'genomes.txt',
               'email': P.PARAMS['hub_email'],
               'descriptionUrl': f'http://userweb.molbiol.ox.ac.uk/{P.PARAMS["hub_publoc"].strip("/")}',
               }

    with open(outfile, 'w') as w:
        for label, info in content.items():
            w.write(f'{label} {info}\n')

@transform(generate_hub_metadata, regex(r'.*.txt'), 'hub_address.txt')
def get_hub_address(infile, outfile):
    with open(outfile, 'w') as w:
        w.write(f'http://userweb.molbiol.ox.ac.uk/{os.path.abspath(infile).lstrip("/")}')

@follows(generate_hub_metadata)
@originate(os.path.join(hub_dir, 'genomes.txt'))
def generate_assembly_metadata(outfile):

    content = {'genome': P.PARAMS['hub_genome'],
               'trackDb': os.path.join(P.PARAMS['hub_genome'], 'trackDb.txt')}

    with open(outfile, 'w') as w:
        for label, info in content.items():
            w.write(f'{label} {info}\n')


@follows(generate_hub_metadata, mkdir(assembly_dir))
@merge([make_bigwig_single, 
        make_bigwig_paired,
        convert_bed_to_bigbed], 
       f'{assembly_dir}/trackDb.txt')
def generate_trackdb_metadata(infiles, outfile):
    def get_track_data(fn):
        return {'track': fn,
                'bigDataUrl': f'http://userweb.molbiol.ox.ac.uk/{(os.path.join(assembly_dir, fn)).lstrip("/")}',
                'shortLabel': fn,
                'longLabel': fn,
                'type': f'{fn.split(".")[-1]}',
                'autoscale': 'on',
                'windowingFunction': 'mean',
                }
     
    # Generate all separate tracks
    bigwig_tracks_all = [get_track_data(os.path.basename(fn)) for fn in infiles]
    
    # Add colours to tracks
    colors = sns.color_palette('husl', len(bigwig_tracks_all))
    for track, color in zip(bigwig_tracks_all, colors):
        track['color'] = ','.join([str(c * 255) for c in color])
    
    
    # Write track data separated
    with open(outfile, 'w') as w:
        for track in bigwig_tracks_all:
            for label, data in track.items():
                w.write(f'{label} {data}\n')
            # Need to separate each track with a new line
            w.write('\n')

@follows(generate_trackdb_metadata)
@transform([make_bigwig_single, make_bigwig_paired, convert_bed_to_bigbed],
           regex(r'(peaks|bigwigs)/(.*)'),
           f'{assembly_dir}/' + r'\2')
def link_hub_files(infile, outfile):
        
    infile_fp = os.path.abspath(infile)
    os.symlink(infile_fp, outfile)



def run_pipeline():
    sys.exit( P.main(sys.argv) )


if __name__ == "__main__":
    sys.exit( P.main(sys.argv) )
    























