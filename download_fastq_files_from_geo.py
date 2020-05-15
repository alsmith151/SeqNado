#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 12:16:20 2019

@author: asmith
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
import requests
import shutil
from bs4 import BeautifulSoup
from Bio import Entrez
import re
from multiprocessing import Pool
import subprocess as sub
import glob
import urllib.request
from datetime import datetime

class NGSRun():

    def __init__(self, accession, label='', force_download=False):
        self.acc = accession
        self.label = label
        self.force = force_download
        self.urls = []
        self.fq_files=[]
        self.metadata = self.get_metadata()
        self.is_paired = False
        
        if len(self.metadata['fastq_ftp']) == 2:
            self.is_paired = True
        

    def _download_from_url(self, url):
        fn_remote = os.path.basename(url)
        if self.label:
            fn_local = fn_remote.replace('.fastq.gz', f'_{self.label}.fastq.gz')
        else:
            fn_local = fn_remote
            
        if not os.path.exists(fn_local) and not self.force:
            print(f'Downloading {fn_remote} to {fn_local}')
            with urllib.request.urlopen(url) as response, open(fn_local, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
        else:
            print(f'File {fn_remote} has already been downloaded -> {fn_local}')
        
        self.fq_files.append(fn_local)
        return fn_local
    
    def download_files(self):
        if self.is_paired:
            for url in self.metadata['fastq_ftp']:
                url_ftp = f'ftp://{url}'
                self._download_from_url(url_ftp)
        else:
            self._download_from_url(f'ftp://{self.metadata["fastq_ftp"][0]}')   
    

    def get_metadata(self):
        
        url = ''.join(['http://www.ebi.ac.uk/ena/data/warehouse/filereport?',
                          f'accession={self.acc}',
                          '&result=read_run&fields=run_accession,fastq_ftp,fastq_md5',
                         ])
    
        metadata = requests.get(url) #Get metadata from EBI ENA warehouse
        
        metadata = metadata.text.strip().split('\n') #Split into column labels and values
        keys, values = [x.split('\t') for x in metadata] 
        metadata = dict(zip(keys, values)) # Convert to dictionary
        metadata = {k:v.split(';') for k, v in metadata.items()} #Separate multiple files
        
        return metadata

class GeoNGS():
    def __init__(self, accession, label, email='asmith151@outlook.com', force=False):
        self.email = email
        self.label = label
        self.force_download = force
        self.acc = accession
        self.run_accs = self.get_srr_acc()
        self.run_objs = [NGSRun(acc, label=self.label, force_download=self.force_download) for acc in self.run_accs]
        
    def get_srr_acc(self):
        
        if not self.email:
            raise ValueError('No email address provided!')
            
        Entrez.email = self.email
        record = Entrez.read(Entrez.esearch('sra', self.acc))
        try:
            uid = record['IdList'][0]
        except Exception as e:
            print(f'Error {e} occured at {self.acc}')
        records = Entrez.esummary(db='sra', id=uid, retmode='xml')
        
        s = BeautifulSoup(records, features='lxml')
        acc_regex = re.compile(r'SRR[0-9]+')
        
        runs = s.find('item', {'name': 'Runs'})
        
        #There appear to be some situations in which multiple runs are present for one GSM number
        # This allows for all of these accessions to be caught.
        return acc_regex.findall(str(runs))
    
    def download_files(self):
        for run in self.run_objs:
            run.download_files()

    def __repr__(self):
        return f'{self.acc}_{self.label}'
    
    def __str__(self):
        return f'{self.acc}_{self.label}'

               
def read_sample_csv(csv):
    df = pd.read_csv(csv, header=None)    
    #Check if csv is comma or tab separated
    if '\t' in ''.join([str(x) for x in df.iloc[0, :]]):
        df = pd.read_csv(csv, sep='\t', header=None)
    
    return (df.values[:, 0], df.values[:, 1])

def download_samples_helper(sample):
    sample.download_files()
    return sample

def generate_md5_checksum_file(samples, fn_md5='md5_checksum_file.txt'):
    '''Uses the metadata associated with each  sample to check the md5 checksums of the downloaded files
       with the supplied md5 checksums'''
    
    
    with open(fn_md5, 'w') as w:       
        for sample in samples:                      # Iterate all of the GEO samples
            for run in sample.run_objs:             # Check each of the associated runs
                if run.is_paired:                   # If the run is paired end then write out the checksums for both files
                    
                    for md5, fn in zip(run.metadata['fastq_md5'], run.fq_files):
                        w.write(f'{md5}  {fn}\n')
                
                else:
                    w.write(f'{run.metadata["fastq_md5"][0]}  {run.fq_files[0]}\n')
        
    
    return fn_md5
                
def compare_md5_checksums(samples, md5_file):
    fn_md5_summary = f'md5_summary_{datetime.now()}.txt'
    cmd = f'md5sum -c {md5_file}'
    sub.run(cmd.split(), stdout=open(fn_md5_summary, 'wb'))
    return pd.read_csv(fn_md5_summary, header=None, sep=' ')



def main():
    p = argparse.ArgumentParser()
    p.add_argument('csv', 
                              help='''Comma or tab delemitated text file. The first column must contain sample accessions (GSM..) 
                              and the second column contains some identifier e.g. assay name''')
    p.add_argument('-p', '--number_of_cores', 
                              help='''Number of cores to use for simultanous sample processing''',
                              default=8, type=int)
    p.add_argument('--force', help='''Forces download of all previously downloaded fq files''',
                              action='store_true')
    
    args = p.parse_args()
    wp = Pool(args.number_of_cores)
    
     
    print('Reading samplesheet')
    samples, labels = read_sample_csv(args.csv)
    print(f'Found {len(samples)} sample(s)')

    print('\nConverting GEO accessions to run accessions stored on EBI ENA')
    geo_samples = [GeoNGS(sample, label) for sample, label in zip(samples, labels)]
    print('Converted sample accessions to run accessions\n')

    for sample in geo_samples:
        print(f'{sample} -> {sample.run_accs}')

    #print(f'\nStarting download for {" , ".join(geo_samples)}')
    geo_samples = wp.map(download_samples_helper, geo_samples)
           
    md5_file = generate_md5_checksum_file(geo_samples)
    md5_df = compare_md5_checksums(geo_samples, md5_file)
    print(md5_df)



if __name__ == '__main__':
    main()

