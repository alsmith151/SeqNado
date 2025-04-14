import datetime
import json
import os
import pathlib
import sys
from typing import Literal, Optional, Dict, List

from jinja2 import Environment, FileSystemLoader
from loguru import logger
from pydantic import BaseModel, field_validator, ValidationError

package_dir = os.path.dirname(os.path.abspath(__file__))
template_dir = os.path.join(package_dir, "workflow/config")

# Define Pydantic Models
class GenomeConfig(BaseModel):
    star_index: str
    bt2_index: str
    chromosome_sizes: str
    gtf: str
    blacklist: Optional[str] = None
    genes: Optional[str] = None
    fasta: Optional[str] = None

class WorkflowConfig(BaseModel):
    # Core Configuration
    assay: Literal["rna", "chip", "atac", "snp", "cat", "meth", 'mcc']
    username: str
    project_date: str
    project_name: str
    seqnado_version: str
    genome: str
    index: str
    chromosome_sizes: str
    gtf: str
    
    # Optional Genome Features
    blacklist: Optional[str] = None
    fastq_screen: bool = False
    fastq_screen_config: Optional[str] = None
    remove_blacklist: bool = False
    
    # PCR Duplicate Handling
    remove_pcr_duplicates: bool = False
    remove_pcr_duplicates_method: Literal["picard", "samtools", "False"] = "False"
    library_complexity: bool = False
    
    # Assay-Specific Features
    shift_atac_reads: bool = False
    spikein: bool = False
    normalisation_method: Optional[Literal["orlando", "with_input"]] = None
    reference_genome: Optional[str] = None
    spikein_genome: Optional[str] = None
    make_bigwigs: bool = False
    pileup_method: Literal["deeptools", "homer", "False"] = "False"
    make_heatmaps: bool = False
    call_peaks: bool = False
    peak_calling_method: Literal["lanceotron", "macs", "homer", "seacr", "False"] = "False"
    consensus_counts: bool = False
    call_methylation: bool = False
    methylation_assay: Literal["bisulfite", "taps", "False"] = "False"
    
    # RNA-Specific
    rna_quantification: Literal["feature_counts", "salmon", "False"] = "False"
    salmon_index: Optional[str] = None
    run_deseq2: bool = False
    
    # SNP-Specific
    call_snps: bool = False
    snp_calling_method: Literal["bcftools", "deepvariant", "False"] = "False"
    snp_database: Optional[str] = None

    # MCC-Specific
    viewpoints: Optional[str] = None
    resolution: Optional[int] = 100

    # For MCC and SNP
    fasta: Optional[str] = None
    fasta_index: Optional[str] = None
    
    # Methylation-Specific
    call_methylation: bool = False
    methylation_assay: Literal["taps", "bisulfite", "False"] = "False"

    # UCSC Hub
    make_ucsc_hub: bool = False
    UCSC_hub_directory: str = "seqnado_output/hub/"
    email: str
    color_by: str = "samplename"
    
    # Additional Options
    options: str  # Raw YAML string for tool options
    geo_submission_files: bool = False
    perform_plotting: bool = False
    plotting_coordinates: Optional[str] = None
    plotting_genes: Optional[str] = None

    @field_validator("remove_pcr_duplicates_method")
    def validate_pcr_method(cls, v, values):
        if values.data.get("remove_pcr_duplicates") and v == "False":
            raise ValueError("Method required when removing duplicates")
        return v

# Helper Functions
def get_user_input(
    prompt: str,
    default: Optional[str] = None,
    is_boolean: bool = False,
    choices: Optional[List[str]] = None,
    is_path: bool = False,
) -> str:
    """
    Prompt the user for input with validation for choices, boolean values, or path existence.
    Re-prompts until valid input is provided.
    """
    while True:
        # Construct the prompt suffix based on choices or default
        if choices:
            prompt_suffix = f"({'/'.join(choices)})"
        elif default is not None:
            prompt_suffix = f"(default: {default})"
        else:
            prompt_suffix = ""
        
        user_input = input(f"{prompt} {prompt_suffix}: ").strip()

        # Handle empty input and apply default if available
        if not user_input:
            if default is not None:
                user_input = default
            else:
                print("Input cannot be empty. Please try again.")
                continue

        # Validate boolean input
        if is_boolean:
            if user_input.lower() in {"yes", "y", "true", "1"}:
                return True
            elif user_input.lower() in {"no", "n", "false", "0"}:
                return False
            else:
                print("Invalid boolean value. Please enter yes/no, y/n, true/false, or 1/0.")
                continue

        # Validate against allowed choices
        if choices and user_input not in choices:
            print(f"Invalid choice. Please choose from: {', '.join(choices)}")
            continue

        # Validate path existence if required
        if is_path and not os.path.exists(user_input):
            print(f"The path '{user_input}' does not exist. Please try again.")
            continue

        return user_input


def load_genome_config() -> Dict[str, GenomeConfig]:
    config_path = pathlib.Path(os.getenv("SEQNADO_CONFIG", pathlib.Path.home())) / ".config/seqnado/genome_config.json"
    if not config_path.exists():
        logger.error("Genome config not found. Run 'seqnado-init' first.")
        sys.exit(1)
    
    with open(config_path) as f:
        return {k: GenomeConfig(**v) for k, v in json.load(f).items()}

def build_workflow_config(assay: str, seqnado_version: str) -> WorkflowConfig:
    genomes = load_genome_config()
    username = os.getenv("USER", "unknown_user")
    today = datetime.datetime.now().strftime("%Y-%m-%d")

    project_name = get_user_input("Project name?", default=f"{username}_project").replace(" ", "_")

    # Validate genome input: list available genomes and re-prompt if needed.
    available_genomes = list(genomes.keys())
    genome = get_user_input(f"Genome? (Available: {', '.join(available_genomes)})", default="hg38")
    while genome not in genomes:
        print(f"Genome '{genome}' is not configured. Please choose from: {', '.join(available_genomes)}")
        genome = get_user_input("Genome?", default="hg38")

    genome_config = genomes[genome]
    index = genome_config.star_index if assay == "rna" else genome_config.bt2_index

    # Base configuration
    config = {
        "assay": assay,
        "username": username,
        "project_date": today,
        "project_name": project_name,
        "seqnado_version": seqnado_version,
        "genome": genome,
        "index": index,
        "chromosome_sizes": genome_config.chromosome_sizes,
        "gtf": genome_config.gtf,
        "blacklist": genome_config.blacklist,
        "email": f"{username}@example.com",
    }

    # Add conditional features
    config.update(get_conditional_features(assay, genome_config))
    try:
        workflow_config = WorkflowConfig(**config)
    except ValidationError as e:
        logger.error(f"Configuration validation error: {e}")
        sys.exit(1)
    return workflow_config

def get_conditional_features(assay: str, genome_config: dict) -> dict:
    features = {}
    username = os.getenv("USER", "unknown_user")

    # Fastq Screen
    features["fastq_screen"] = get_user_input("Perform fastqscreen?", default="no", is_boolean=True)
    if features["fastq_screen"]:
        features["fastq_screen_config"] = get_user_input("Fastqscreen config path:", default="/ceph/project/milne_group/shared/seqnado_reference/fastqscreen_reference/fastq_screen.conf", is_path=True)
    
    # Blacklist Handling
    features["remove_blacklist"] = get_user_input("Remove blacklist regions?", default="yes", is_boolean=True)
    
    # PCR Duplicates
    default_duplicates = "yes" if assay in ["chip", "atac", "cat", "meth"] else "no"
    features["remove_pcr_duplicates"] = get_user_input("Remove PCR duplicates?", default=default_duplicates, is_boolean=True)
    if features["remove_pcr_duplicates"]:
        features["remove_pcr_duplicates_method"] = get_user_input("Duplicates removal method:", choices=["picard", "samtools"], default="picard")
        features["library_complexity"] = get_user_input("Calculate library complexity?", default="no", is_boolean=True)
    
    # ATAC-Specific Logic
    if assay == "atac":
        features["shift_atac_reads"] = get_user_input("Shift ATAC reads?", default="yes", is_boolean=True)
    
    # MCC-Specific Logic
    if assay == "mcc":
        features['viewpoints'] = get_user_input("Path to viewpoints file:", default="path/to/viewpoints.bed", is_path=False)
        features['fasta'] = get_user_input("Path to reference fasta:", default="path/to/reference.fasta", is_path=False)
        features['resolution'] = get_user_input("Resolution for MCC cooler files:", default="100")
    
    # Spike-in Normalisation
    if assay in ["chip", "rna", 'cat']:
        features["spikein"] = get_user_input(
            "Do you have spikein? (yes/no)", default="no", is_boolean=True
        )
        if features["spikein"] and not assay == "rna":
            features["normalisation_method"] = get_user_input(
                "Normalisation method:",
                default="orlando",
                choices=["orlando", "with_input"],
            )
            features["reference_genome"] = get_user_input(
                "Reference genome:", default="hg38"
            )
            features["spikein_genome"] = get_user_input(
                "Spikein genome:", default="dm6"
            )
    
    # RNA-Specific Features
    if assay == "rna":
        features["rna_quantification"] = get_user_input("Quantification method:", choices=["feature_counts", "salmon"], default="feature_counts")
        if features["rna_quantification"] == "salmon":
            features["salmon_index"] = get_user_input("Salmon index path:", default="path/to/salmon_index")
        features["run_deseq2"] = get_user_input("Run DESeq2?", default="no", is_boolean=True)
    
    # Peak calling options
    if assay in ["chip", "atac", "cat"]:
        features["call_peaks"] = get_user_input("Call peaks?", default="yes", is_boolean=True)
        if features["call_peaks"]:
            
            match assay:
                case "chip":
                    default = "lanceotron"
                case "atac":
                    default = "lanceotron"
                case "cat":
                    default = 'seacr'
                case _:
                    default = "lanceotron"

            features["peak_calling_method"] = get_user_input("Peak calling method:", choices=["lanceotron", "macs", "homer", "seacr"], default=default)
    
    # Pileup method
    if assay != "snp":
        features['make_bigwigs'] = get_user_input("Make Bigwigs?", default="no", is_boolean=True)
        if features['make_bigwigs']:
            features['pileup_method'] = get_user_input("Bigwig method:", choices=["deeptools", "homer"], default="deeptools")
    
    # Heatmaps
    features["make_heatmaps"] = get_user_input("Make heatmaps?", default="no", is_boolean=True)
    
    # UCSC Hub
    features["make_ucsc_hub"] = get_user_input("Make UCSC hub?", default="no", is_boolean=True)
    if features["make_ucsc_hub"]:
        features["UCSC_hub_directory"] = get_user_input("UCSC hub directory:", default="seqnado_output/hub/")
        features["email"] = get_user_input("What is your email address?", default=f"{username}@example.com")
        features["color_by"] = get_user_input("Color by (for UCSC hub):", default="samplename")
    
    # SNP Calling
    if assay == "snp":
        features["call_snps"] = get_user_input("Call SNPs?", default="no", is_boolean=True)
        if features["call_snps"]:
            features["snp_calling_method"] = get_user_input("SNP caller:", choices=["bcftools", "deepvariant"], default="bcftools")
            features["fasta"] = genome_config.fasta if genome_config.fasta else get_user_input("Path to reference fasta:", default='no')
            features["fasta_index"] = get_user_input("Path to reference fasta index:", default="path/to/reference.fasta.fai", is_path=True)
            features["snp_database"] = get_user_input("Path to SNP database:", default="path/to/snp_database", is_path=True)
    
    # Methylation Calling
    if assay == "meth":
        features["call_methylation"] = get_user_input("Call methylation?", default="no", is_boolean=True)
        if features["call_methylation"]:
            features["fasta"] = genome_config.fasta if genome_config.fasta else get_user_input("Path to reference fasta:", default='no')
            features["methylation_assay"] = get_user_input("Methylation assay:", choices=["taps", "bisulfite"], default="taps")

    # Consensus Counts
    if not assay == 'rna':
        features["consensus_counts"] = get_user_input(
            "Generate consensus counts from Design merge column? (yes/no)",
            default="no",
            is_boolean=True,
        )

    # GEO Submission
    features["geo_submission_files"] = get_user_input("Generate GEO submission files?", default="no", is_boolean=True)
    
    # Plotting
    features["perform_plotting"] = get_user_input("Perform plotting?", default="no", is_boolean=True)
    if features["perform_plotting"]:
        features["plotting_coordinates"] = get_user_input("Path to bed file with coordinates for plotting", default=None)
        features["plotting_genes"] = genome_config.genes if genome_config.genes else get_user_input("Path to bed file with genes.", default='no')
    
    # Add tool options
    features["options"] = get_tool_options(assay)
    return features

def get_tool_options(assay: str) -> str:
    """
    Return the tool options YAML string for the given assay.

    Args:
        assay (str): The assay type, used to determine the correct tool options.
    
    Returns:
        str: The YAML string with tool options for the given
    """

    import importlib.resources
    import yaml
    import seqnado.workflow.config

    match assay:
        case "rna":
            tool_file = importlib.resources.files(seqnado.workflow.config) / 'tool_options_rna.yml'
        case "snp":
            tool_file = importlib.resources.files(seqnado.workflow.config) / 'tool_options_snp.yml'
        case "meth":
            tool_file = importlib.resources.files(seqnado.workflow.config) / 'tool_options_meth.yml'
        case "mcc":
            tool_file = importlib.resources.files(seqnado.workflow.config) / 'tool_options_mcc.yml'
        case _:
            tool_file = importlib.resources.files(seqnado.workflow.config) / 'tool_options_base.yml'
    

    with open(tool_file) as f:
        tool_options = yaml.safe_load(f)


    # Make assay specific changes
    if assay == "mcc":
        tool_options['samtools']['filter_options'] = ''
        tool_options['deeptools']['bamcoverage'] = ''


    return yaml.dump(tool_options)



# Template Rendering
def create_config(assay: str, rerun: bool, seqnado_version: str, debug=False):
    env = Environment(loader=FileSystemLoader(template_dir))
    workflow_config = build_workflow_config(assay, seqnado_version)
    
    dir_name = os.getcwd() if rerun else f"{workflow_config.project_date}_{workflow_config.assay}_{workflow_config.project_name}"
    os.makedirs(dir_name, exist_ok=True)

    fastq_dir = pathlib.Path(dir_name) / "fastq"
    fastq_dir.mkdir(exist_ok=True)

    # Render main config
    with open(f"{dir_name}/config_{assay}.yml", "w") as f:
        f.write(env.get_template("config.yaml.jinja").render(workflow_config.model_dump()))
    
    # Additional RNA template
    if assay == "rna":
        with open(f"{dir_name}/deseq2_{workflow_config.project_name}.qmd", "w") as f:
            f.write(env.get_template("deseq2.qmd.jinja").render(workflow_config.model_dump()))

    logger.success(f"Created configuration in {dir_name}")

# # Preserve original tool option YAML strings
# TOOL_OPTIONS = """
# trim_galore:
#     threads: 4
#     options: --2colour 20 

# bowtie2:
#     threads: 8
#     options:

# samtools:
#     threads: 16
#     filter_options: -f 2

# picard:
#     threads: 8
#     options:

# homer:
#     use_input: true
#     maketagdirectory:
#     makebigwig:
#     findpeaks:

# deeptools:
#     threads: 8
#     alignmentsieve: --minMappingQuality 30 
#     bamcoverage: --extendReads -bs 1 --normalizeUsing RPKM --minMappingQuality 10

# macs:
#     version: 2
#     callpeak: -f BAMPE

# lanceotron:
#     use_input: True
#     callpeak: -c 0.5

# seacr:
#     threshold: 0.01
#     norm: non
#     stringency: stringent

# heatmap:
#     options: -b 1000 -m 5000 -a 1000 --binSize 50
#     colormap: RdYlBu_r 

# featurecounts:
#     threads: 16
#     options:  -p --countReadPairs
    
# """

# TOOL_OPTIONS_RNA = """
# trim_galore:
#     threads: 4
#     options: --2colour 20 

# star:
#     threads: 16
#     options: --quantMode TranscriptomeSAM GeneCounts --outSAMunmapped Within --outSAMattributes Standard

# samtools:
#     threads: 16
#     filter_options: -f 2

# picard:
#     threads: 8
#     options:

# featurecounts:
#     threads: 16
#     options: -s 0 -p --countReadPairs -t exon -g gene_id

# salmon:
#     threads: 16
#     options: --libType A
    
# homer:
#     maketagdirectory:
#     makebigwig:

# deeptools:
#     threads: 16
#     alignmentsieve: --minMappingQuality 30 
#     bamcoverage: -bs 1 --normalizeUsing CPM

# heatmap:
#     options: -b 1000 -m 5000 -a 1000 --binSize 50
#     colormap: RdYlBu_r 
# """

# TOOL_OPTIONS_SNP = """
# trim_galore:
#     threads: 8
#     options: --2colour 20 

# bowtie2:
#     threads: 8
#     options:

# samtools:
#     threads: 16
#     filter_options: -f 2
    
# picard:
#     threads: 8
#     options:

# bcftools:
#     threads: 16
#     options:
    
# """

# TOOL_OPTIONS_METH = """
# trim_galore:
#     threads: 8
#     options: --2colour 20 

# bowtie2:
#     threads: 8
#     options:

# samtools:
#     threads: 16
#     filter_options: -f 2
    
# picard:
#     threads: 8
#     options:

# methyldackel:
#     threads: 16
#     options:
    
# """