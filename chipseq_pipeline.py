#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 10:47:03 2020

@author: asmith
"""

# import packages
import sys
import os
import gzip
import seaborn as sns
from cgatcore import pipeline as P
from ruffus import mkdir, follows, transform, merge, originate, collate, split, regex, add_inputs, suffix

# Read in parameter file
P.get_parameters('chipseq_pipeline.yml')

# Global variables
hub_dir = os.path.join(P.PARAMS["hub_publoc"], P.PARAMS['hub_name'])
assembly_dir = os.path.join(hub_dir, P.PARAMS['hub_genome'])


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
          job_threads=P.PARAMS['threads'])


@follows(mkdir('report'))
@merge(qc_reads, 'report/readqc_report.html')
def multiqc_reads (infile, outfile):
    '''Collate fastqc reports into single report using multiqc'''
    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc fastqc/ -o report -n readqc_report.html'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_memory='16G')
    
    
@follows(mkdir('trimmed'))
@transform('*.fastq.gz', 
           regex(r'(?!.*_.*_[12])^(.*).fastq.gz'), # Regex negates any filenames matching the paired pattern
           r'trimmed/\1_trimmed.fq.gz')
def trim_reads_single(infile, outfile):
    
    statement = '''trim_galore --cores %(threads)s %(trim_options)s 
                  -o trimmed %(infile)s'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'])

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
          job_threads=P.PARAMS['threads'])


@follows(mkdir('bam'), trim_reads_single)
@transform(trim_reads_single, 
           regex(r'trimmed/(.*)_trimmed.fastq.gz'), 
           r'bam/\1.bam')
def align_reads_single(infile, outfile):
    ''' Aligns digested fq files using bowtie2'''
    
    options = P.PARAMS['bowtie2_options'] if P.PARAMS['bowtie2_options'] else ''
        
    statement = '''bowtie2 -x %(bowtie2_index)s -U %(infile)s 
                    -p %(threads)s %(options)s 
                    | samtools view -bS > %(outfile)s 2> %(outfile)s.log
                    && samtools sort %(outfile)s -o %(outfile)s.sorted.bam -m 2G -@ %(threads)s
                    && mv %(outfile)s.sorted.bam %(outfile)s'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'],
          job_memory='20G')

   
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
        
    if P.PARAMS['bowtie2_options']:
        options = P.PARAMS['bowtie2_options']
    
    cmd = '''bowtie2 -x %(bowtie2_index)s -1 %(fq1)s -2 %(fq2)s -p %(threads)s %(options)s |
             samtools view -b - > %(outfile)s &&
             samtools sort -@ %(threads)s -m 5G -o %(sorted_bam)s %(outfile)s &&
             mv %(sorted_bam)s %(outfile)s'''

    P.run(cmd, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'])


@transform([align_reads_single, align_reads_paired], 
           regex(r'bam/(.*).bam'),
           r'bam/\1.bam.bai')
def create_bam_file_index(infile, outfile):
    """Bam files are compressed. The index allows fast access to different
    slices of the file."""
    statement = 'samtools index %(infile)s %(outfile)s'
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'])


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
          job_queue=P.PARAMS['queue'])

@merge(mapping_qc, 'report/mapping_report.html')
def mapping_multiqc (infile, outfile):
    '''Combines mapping metrics using multiqc'''
    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc bam/ -o report -n mapping_report.html'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_memory='16G')


@follows(create_bam_file_index, mkdir('deduplicated'))    
@transform([align_reads_single, align_reads_paired],
           regex(r'bam/(.*.bam)'),
           r'deduplicated/\1')
def remove_duplicates(infile, outfile):
    stats = outfile.replace('.bam', '_marked_duplicates.txt')
    statement = ' '.join(['picard MarkDuplicates',
                          'I=%(infile)s',
                          'O=%(outfile)s',
                          'M=%(stats)s',
                          'REMOVE_DUPLICATES=true',
                          'REMOVE_SEQUENCING_DUPLICATES=true',
                          'CREATE_INDEX=true'])
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'])



@follows(mkdir('bigwigs'), create_bam_file_index, remove_duplicates)
@transform(align_reads_paired, regex(r'bam/(.*).bam'), r'bigwigs/\1.bigWig')
def make_bigwig_paired(infile, outfile):
    
    cmd = [ 'bamCoverage', 
            '-b', infile, 
            '-o', outfile,
            '-p', '%(threads)s',
            '--effectiveGenomeSize', '%(genome_size)s',
            '--normalizeUsing', 'BPM',
            '--smoothLength', '%(bigwig_smoothing_window)s',
            '--extendReads',
            '--verbose']
    
    statement = ' '.join(cmd)
    
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'],
          job_threads=P.PARAMS['threads'],)
    
@follows(mkdir('bigwigs'), create_bam_file_index, remove_duplicates)
@transform(align_reads_single, regex(r'bam/(.*).bam'), r'bigwigs/\1.bigWig')
def make_bigwig_single(infile, outfile):
    
    cmd = [ 'bamCoverage', 
            '-b', infile, 
            '-o', outfile,
            '-p', '%(threads)s',
            '--effectiveGenomeSize', '%(genome_size)s',
            '--normalizeUsing', 'BPM',
            '--smoothLength', '%(bigwig_smoothing_window)s',
            '--verbose']
    
    statement = ' '.join(cmd)
    
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'],
          job_threads=P.PARAMS['threads'],)   
    


@follows(mkdir('peaks'))
@transform(remove_duplicates, 
           regex(r'deduplicated/(.*)_(?!input|Input|INPUT)(.*).bam'),
           add_inputs(r'deduplicated/\1_input.bam'),
           r'peaks/\1_\2_peaks.narrowPeak')
def call_peaks(infile, outfile):
    
    treatment = infile[0]
    control= infile[1]
    
    file_base = outfile.replace('_peaks.narrowPeak', '')
    macs2_options = P.PARAMS['macs2_options'] if P.PARAMS['macs2_options'] else ''
    
    statement = '''macs2 callpeak %(macs2_options)s -g %(genome_size)s 
                  -t %(treatment)s -n %(file_base)s'''
    
    if control:
        statement += ' -c %(control)s'
      

    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'])

@transform(call_peaks, regex(r'peaks/(.*).narrowPeak'), r'peaks/\1.bed')
def convert_narrowpeak_to_bed(infile, outfile):
    
    statement = '''awk '{OFS="\\t"; print $1,$2,$3,$4}' %(infile)s > %(outfile)s'''
    
    P.run(statement,
          job_queue  = P.PARAMS['queue'])
    
@transform(convert_narrowpeak_to_bed, regex(r'peaks/(.*).bed'), r'peaks/\1.bigBed')    
def convert_bed_to_bigbed(infile, outfile):
    
    statement = '''bedToBigBed %(infile)s %(genome_chrom_sizes)s %(outfile)s'''
    P.run(statement,
          job_queue  = P.PARAMS['queue'])
        
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


if __name__ == "__main__":
    sys.exit( P.main(sys.argv) )























