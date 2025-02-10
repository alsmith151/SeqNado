import pandas as pd
import subprocess
import re

assay = snakemake.params.assay
config = snakemake.config


def fastqc_version():
    cmd = 'fastqc --version'
    version = subprocess.check_output(cmd, shell=True).decode().strip()
    version = re.search(r'v(\d+).(\d+).(\d+)', version).group(0)
    return version

def trim_galore_version():
    cmd = 'trim_galore --version'
    version = subprocess.check_output(cmd, shell=True).decode().strip()
    version = re.search(r'version (\d+).(\d+).(\d+)', version).group(0)
    return version

def star_version():
    cmd = 'STAR --version'
    version = subprocess.check_output(cmd, shell=True).decode().strip()
    return version

def bowtie2_version():
    cmd = 'bowtie2 --version'
    version = subprocess.check_output(cmd, shell=True).decode().strip()
    version = re.search(r'bowtie2-align-s version (\d+).(\d+).(\d+)', version).group(0)
    return version

def featureCounts_version():
    cmd = 'featureCounts -v'
    version = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode().strip().replace('featureCounts v', '')
    return version

def bamCoverage_version():
    cmd = 'bamCoverage --version'
    version = subprocess.check_output(cmd, shell=True).decode().strip()
    return version

def picard_version():
    cmd = 'picard MarkDuplicates --version'
    try:
        version = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode().strip().removeprefix('Version:')
    except subprocess.CalledProcessError:
        version = "3.1.1"
    return version

def lanceotron_version():
    cmd = 'lanceotron --version'
    version = subprocess.check_output(cmd, shell=True).decode().strip()
    return version

def macs2_version():
    cmd = 'macs2 --version'
    version = subprocess.check_output(cmd, shell=True).decode().strip()
    re.search(r'(\d|\.)+$', version).group(1)
    return version




content = f"""
FASTQ files were quality checked using FastQC version {fastqc_version()}.
Adapter sequences were removed and low quality reads were trimmed using trim_galore {trim_galore_version()} using the following parameters: {config['trim_galore']['options']}.
Reads were aligned to the reference genome {config['genome']['name']} using {'STAR' if assay == 'RNA' else 'bowtie2'} v{star_version() if assay == 'RNA' else bowtie2_version()} with the following parameters: {config['star' if assay == 'RNA' else 'bowtie2']['options']}.
"""

if config['remove_pcr_duplicates_method']:
    content += f"""Duplicate reads were removed using Picard MarkDuplicates v{picard_version()} with the following parameters: {config['picard']['options']}"""


content += f"""
BigWig files were generated using deepTools {bamCoverage_version()} with the following parameters: {config['deeptools']['bamcoverage']}.
"""

if assay == 'RNA':
    content += f"""Strands were separated using the --filterRNAstrand option. """
    content += f"""Alignments were quantified using featureCounts v{featureCounts_version()} with the following parameters: {config['featurecounts']['options']}"""


if assay in ['ChIP', 'ATAC', 'CUT&TAG']:
    if "lanceotron" in  config['peak_calling_method']:
        content += f"""Peak calling was performed using lanceotron v1.2.6 with the following parameters: {config['lanceotron']['callpeak']}"""
    
    if 'macs2' in config['peak_calling_method']:
        content += f"""Peak calling was performed using MACS2 v{macs2_version()} with the following parameters: {config['macs2']['callpeak']}"""



content = content.strip().replace('the following parameters: False', 'with default parameters').replace('/n/n', '/n')
content = "\n".join([l.strip('\n') for l in content.splitlines() if l.strip()])


with open(snakemake.output[0], 'w') as f:
    f.write(content)

