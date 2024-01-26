import os
import datetime
from jinja2 import Environment, FileSystemLoader
import json

package_dir = os.path.dirname(os.path.abspath(__file__))
template_dir = os.path.join(package_dir, "workflow/config")

def get_user_input(prompt, default=None, is_boolean=False, choices=None):
    while True:
        user_input = input(f"{prompt} [{'/'.join(choices) if choices else default}]: ") or default
        if is_boolean:
            return user_input.lower() == 'yes'
        if choices and user_input not in choices:
            print(f"Invalid choice. Please choose from {', '.join(choices)}.")
            continue
        return user_input

def load_genome_values():
    with open(os.path.join(template_dir, 'preset_genomes.json'), 'r') as f:
        return json.load(f)

def prompt_for_genome_info(assay, genome):
    prompts = {
        "indices": "Path to Bowtie2 genome index:" if assay in ["chip", "atac"] else "Path to STAR v2.7.10b genome index:",
        "chromosome_sizes": "Path to chromosome sizes file:",
        "gtf": "Path to GTF file:",
        "blacklist": "Path to blacklist bed file:"
    }
    return {key: get_user_input(prompt) for key, prompt in prompts.items()}

def get_genome_dict(assay, genome, genome_values):
    if genome == "other":
        genome = get_user_input("What is your genome name?", default="other")
        return {genome: prompt_for_genome_info(assay, genome)}
    elif genome in genome_values:
        index_key = 'bt2_index' if assay in ["chip", "atac"] else 'star_index'
        genome_info = genome_values[genome].copy()
        genome_info['indices'] = genome_info.pop(index_key)
        return {genome: genome_info}
    return {}

def update_template_data(assay, genome_dict, genome, template_data):
    genome_info = genome_dict.get(genome, {})
    
    for key in ['indices', 'chromosome_sizes', 'gtf', 'blacklist']:
        template_data[key] = genome_info.get(key)

    # Common settings for all assays
    template_data['remove_blacklist'] = get_user_input("Do you want to remove blacklist regions? (yes/no)", default="yes", is_boolean=True)
    if template_data['remove_blacklist']:
        template_data['blacklist'] = genome_info.get('blacklist')

    # Assay-specific configurations
    if assay in ["chip", "atac"]:
        template_data['remove_pcr_duplicates'] = get_user_input("Remove PCR duplicates? (yes/no)", default="yes", is_boolean=True)
        
        template_data['call_peaks'] = get_user_input("Do you want to call peaks? (yes/no)", default="no", is_boolean=True)

        if template_data['call_peaks']:
            template_data['peak_calling_method'] = get_user_input("Peak caller:", default="macs", choices=["lanceotron", "macs", "homer"])

        if assay == "chip":
            template_data['spikein'] = get_user_input("Do you have chip spikein? (yes/no)", default="no", is_boolean=True)
            if template_data['spikein']:
                template_data['normalisation_method'] = get_user_input("Normalisation method:", default="orlando", choices=["orlando", "with_input"])
                template_data['reference_genome'] = get_user_input("Reference genome:", default="hg38")
                template_data['spikein_genome'] = get_user_input("Spikein genome:", default="dm6")
                
        if assay == "atac":
            template_data['shift_atac_reads'] = get_user_input("Shift ATAC-seq reads? (yes/no)", default="yes" if assay == "atac" else "no", is_boolean=True)

    elif assay == "rna":
        template_data['remove_pcr_duplicates'] = get_user_input("Remove PCR duplicates? (yes/no)", default="no", is_boolean=True)
        template_data['run_deseq2'] = get_user_input("Run DESeq2? (yes/no)", default="no", is_boolean=True)

    # PCR duplicates method
    if template_data.get('remove_pcr_duplicates'):
        template_data['remove_pcr_duplicates_method'] = get_user_input("Remove PCR duplicates method:", default="picard", choices=["picard"])
    else:
        template_data['remove_pcr_duplicates_method'] = "False"

    # General options for all assays
    template_data['make_bigwigs'] = get_user_input("Do you want to make bigwigs? (yes/no)", default="no", is_boolean=True)
    if template_data['make_bigwigs']:
        template_data['pileup_method'] = get_user_input("Pileup method:", default="deeptools", choices=["deeptools", "homer"])

    template_data['make_ucsc_hub'] = get_user_input("Do you want to make a UCSC hub? (yes/no)", default="no", is_boolean=True)
    if template_data['make_ucsc_hub']:
        template_data['UCSC_hub_directory'] = get_user_input("UCSC hub directory:", default="/path/to/ucsc_hub/")
        template_data['email'] = get_user_input("What is your email address?", default=f"{template_data['username']}@example.com")
        template_data['color_by'] = get_user_input("Color by (for UCSC hub):", default="samplename")
    else:
        template_data['UCSC_hub_directory'] = "."
        template_data['email'] = f"{template_data['username']}@example.com"
        template_data['color_by'] = "samplename"
        
    if assay in ["chip", "atac"]:
        template_data['options'] = TOOL_OPTIONS
    elif assay == "rna":
        template_data['options'] = TOOL_OPTIONS_RNA


def setup_configuration(assay, genome, template_data):
    # Collect common configurations
    username = os.getenv('USER', 'unknown_user')
    today = datetime.datetime.now().strftime('%Y-%m-%d')
    project_name = get_user_input("What is your project name?", default=f"{username}_project")

    common_config = {
        'username': username,
        'project_date': today,
        'project_name': project_name,
        'genome': genome
    }
    template_data.update(common_config)

    # Load genome values
    genome_values = load_genome_values()
    
    # Get genome-specific information
    genome_dict = get_genome_dict(assay, genome, genome_values)
    
    # Update template data with genome-specific and assay-specific information
    update_template_data(assay, genome_dict, genome, template_data)


# Tool Specific Options
TOOL_OPTIONS = """
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
    bamcoverage: --extendReads -bs 1 --normalizeUsing RPKM

macs:
    version: 2
    callpeak:

lanceotron:
    use_input: True
    callpeak: -c 0.5

heatmap:
    options: -m 10000 -b 3000 -a 3000
    colormap: RdYlBu_r
"""

TOOL_OPTIONS_RNA = """
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
    options: -s 0 -p --countReadPairs -t exon -g gene_id

homer:
    maketagdirectory:
    makebigwig:

deeptools:
    threads: 8
    alignmentsieve: --minMappingQuality 30 
    bamcoverage: -bs 1 --normalizeUsing CPM

heatmap:
    options: -m 10000 -b 3000 -a 3000
    colormap: RdYlBu_r
"""

def create_config(assay, genome):
    env = Environment(loader=FileSystemLoader(template_dir), auto_reload=False)
    template = env.get_template("config.yaml.jinja")        
    template_deseq2 = env.get_template("deseq2.qmd.jinja")

    template_data = {'assay': assay, 'genome': genome}
    setup_configuration(assay, genome, template_data)

    dir_name = f"{template_data['project_date']}_{template_data['assay']}_{template_data['project_name']}"
    os.makedirs(dir_name, exist_ok=True)
    fastq_dir = os.path.join(dir_name, "fastq")
    os.makedirs(fastq_dir, exist_ok=True)

    with open(os.path.join(dir_name, f"config_{assay}.yml"), 'w') as file:
        file.write(template.render(template_data))

    if assay == "rna":
        with open(os.path.join(dir_name, f"deseq2_{template_data['project_name']}.qmd"), 'w') as file:
            file.write(template_deseq2.render(template_data))

    print(f"Directory '{dir_name}' has been created with the 'config_{assay}.yml' file.")
