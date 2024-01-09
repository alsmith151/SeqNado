import os
import datetime
from jinja2 import Environment, FileSystemLoader

# Helper Functions
def get_user_input(prompt, default=None, is_boolean=False, choices=None):
    while True:
        user_input = input(f"{prompt} [{'/'.join(choices) if choices else default}]: ") or default
        if is_boolean:
            return user_input.lower() == 'yes'
        if choices and user_input not in choices:
            print(f"Invalid choice. Please choose from {', '.join(choices)}.")
            continue
        return user_input

def setup_configuration(assay, template_data):
    username = os.getenv('USER', 'unknown_user')
    today = datetime.datetime.now().strftime('%Y-%m-%d')
    project_name = get_user_input("What is your project name?", default=f"{username}_project")
    genome = get_user_input("What is your genome name?", default="hg38", choices=['dm6', 'hg19', 'hg38', 'hg38_dm6', 'hg38_mm39', 'hg38_spikein', 'mm9', 'mm10', 'mm10_TetR_ChrX_and_Chr8', 'mm39', 'other'])
    common_config = {
        'username': username,
        'project_date': today,
        'project_name': project_name,
        'genome': genome
    }
    
    template_data.update(common_config)
    
    if assay in ["chip", "atac"]:
        template_data['indicies'] = get_user_input("Path to Bowtie2 genome indices:", default=f"/ceph/project/milne_group/shared/seqnado_reference/{genome}/UCSC/bt2_index/{genome}") 
    elif assay == "rna":
        template_data['indicies'] = get_user_input("Path to STAR genome indices:", default=f"/ceph/project/milne_group/shared/seqnado_reference/{genome}/UCSC/STAR_2.7.10b")

    template_data['chromosome_sizes'] = get_user_input("Path to chromosome sizes file:", default=f"/ceph/project/milne_group/shared/seqnado_reference/{genome}/UCSC/sequence/{genome}.chrom.sizes")
    template_data['gtf'] = get_user_input("Path to GTF file:", default=f"/ceph/project/milne_group/shared/seqnado_reference/{genome}/UCSC/genes/{genome}.ncbiRefSeq.gtf")
    template_data['read_type'] = get_user_input("What is your read type?", default="paired", choices=["paired", "single"])

    template_data['split_fastq'] = get_user_input("Do you want to split FASTQ files? (yes/no)", default="no", is_boolean=True)
    if template_data['split_fastq']:
        template_data.update['split_fastq_parts'] = get_user_input("How many parts do you want to split the FASTQ files into?", default="4")

    if assay in ["chip", "atac"]:
        template_data['remove_pcr_duplicates'] = get_user_input("Remove PCR duplicates? (yes/no)", default="yes", is_boolean=True),
    elif assay == "rna":
        template_data['remove_pcr_duplicates'] = get_user_input("Remove PCR duplicates? (yes/no)", default="no", is_boolean=True),
    
    if template_data['remove_pcr_duplicates']:
        template_data.update({
            'remove_pcr_duplicates_method': get_user_input("Remove PCR duplicates method:", default="picard", choices=["picard", "deeptools"])
        })

    template_data['remove_blacklist'] = get_user_input("Do you want to remove blacklist regions? (yes/no)", default="yes", is_boolean=True)
    if template_data['remove_blacklist']:
        template_data['blacklist'] = get_user_input("Path to blacklist file:", default=f"/ceph/project/milne_group/shared/seqnado_reference/{genome}/UCSC/blacklist/{genome}.blacklist.bed")
    
    if assay == "atac":
        template_data['shift_atac_reads'] = get_user_input("Shift ATAC-seq reads? (yes/no)", default="yes", is_boolean=True)
    elif assay in ["chip", "rna"]:
        template_data['shift_atac_reads'] = "False"

    if assay == "chip":
        template_data['chip_spikein_normalisation'] = get_user_input("Do you have spikein? (yes/no)", default="no", is_boolean=True)
    elif assay in ["atac", "rna"]:
        template_data['chip_spikein_normalisation'] = "False"

    if assay == "chip":
        if template_data['chip_spikein_normalisation']:
                template_data['reference_genome'] = get_user_input("Reference genome:", default="hg38")
                template_data['spikein_genome'] = get_user_input("Spikein genome:", default="dm6")
                template_data['fastq_screen_config'] = get_user_input("Path to fastqscreen config:", default="/ceph/project/milne_group/shared/seqnado_reference/fastqscreen_reference/fastq_screen.conf")
    
    template_data['make_bigwigs'] = get_user_input("Do you want to make bigwigs? (yes/no)", default="yes", is_boolean=True)
    if template_data['make_bigwigs']:
        template_data['pileup_method'] = get_user_input("Pileup method:", default="deeptools", choices=["deeptools", "homer"])
        template_data['make_heatmaps'] = get_user_input("Do you want to make heatmaps? (yes/no)", default="yes", is_boolean=True)
    
    if assay in ["chip", "atac"]:
        template_data['call_peaks'] = get_user_input("Do you want to call peaks? (yes/no)", default="yes", is_boolean=True)
    elif assay == "rna":
        template_data['call_peaks'] = "False"
    if assay in ["chip", "atac"]:
        if template_data['call_peaks']:
            template_data['peak_calling_method'] = get_user_input("Peak caller:", default="macs", choices=["macs", "homer", "lanceotron"])

    if assay == "rna":
        template_data['run_deseq2'] = get_user_input("Run DESeq2? (yes/no)", default="yes", is_boolean=True)
    elif assay in ["chip", "atac"]:
        template_data['run_deseq2'] = "False"

    template_data['make_ucsc_hub'] = get_user_input("Do you want to make a UCSC hub? (yes/no)", default="yes", is_boolean=True)
    if template_data['make_ucsc_hub']:
        template_data['UCSC_hub_directory'] = get_user_input("UCSC hub directory:", default="/path/to/ucsc_hub/")
        template_data['email'] = get_user_input("What is your email address?", default=f"{username}@example.com")
        template_data['color_by'] = get_user_input("Color by (for UCSC hub):", default="samplename")

    if assay in ["chip", "atac"]:
        template_data['options'] = TOOL_OPTIONS
    elif assay == "rna":
        template_data['options'] = TOOL_OPTIONS_RNA


# Tool Specific Options
TOOL_OPTIONS = """
#################################
# Tool specific options         #
#################################

trim_galore:
threads: 4
options: --2colour 20 

bowtie2:
threads: 4
options:

picard:
threads: 4
options:

homer:
use_input: true
maketagdirectory:
makebigwig:
findpeaks:

deeptools:
threads: 8
alignmentsieve: --minMappingQuality 30 
# Options passed to deeptools BamCoverage
# These need to be replaced
# e.g. --extendReads -bs 1 --normalizeUsing RPKM
bamcoverage: --extendReads -bs 1 --normalizeUsing RPKM

macs:
version: 2
callpeak:

lanceotron:
# Instructs lanceotron to use the matched input file for peak calling
# No effect if input file is not matched
use_input: True
# Options passed to callPeaks[Input] command
callpeak: -c 0.5

heatmap:
options:
colormap: RdYlBu_r
"""

TOOL_OPTIONS_RNA = """
#################################
# Tool specific options         #
#################################

trim_galore:
threads: 4
options: --2colour 20 

star:
threads: 4
options: --quantMode TranscriptomeSAM GeneCounts --outSAMunmapped Within --outSAMattributes Standard

picard:
threads: 4
options:

featurecounts:
threads: 4
options: -p

homer:
maketagdirectory:
makebigwig:

deeptools:
threads: 8
alignmentsieve: --minMappingQuality 30 
# Options passed to deeptools BamCoverage
# These need to be replaced
# e.g. --extendReads -bs 1 --normalizeUsing RPKM
bamcoverage: -bs 1 --normalizeUsing CPM
"""

def create_config(assay):
    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template("seqnado/workflow/config/config.yaml.jinja")
    template_deseq2 = env.get_template("seqnado/workflow/config/deseq2.qmd.jinja")
    
    # Initialize template data
    template_data = {'assay': assay}

    # Setup configuration
    setup_configuration(assay, template_data)

    # Create directory and render template
    dir_name = f"{template_data['project_date']}_{template_data['assay']}_{template_data['project_name']}"
    os.makedirs(dir_name, exist_ok=True)
    fastq_dir = os.path.join(dir_name, "fastq")
    os.makedirs(fastq_dir, exist_ok=True)
    
    with open(os.path.join(dir_name, f"config_{assay}.yml"), 'w') as file:
        file.write(template.render(template_data))

    # add deseq2 qmd file if rna
    if assay == "rna":
        with open(os.path.join(dir_name, f"deseq2_{template_data['project_name']}.qmd"), 'w') as file:
            file.write(template_deseq2.render(template_data))
            
    print(f"Directory '{dir_name}' has been created with the 'config_{assay}.yml' file.")

