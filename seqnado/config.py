import os
import datetime
from jinja2 import Environment, FileSystemLoader
import json

package_dir = os.path.dirname(os.path.abspath(__file__))
template_dir = os.path.join(package_dir, "workflow/config")


# Helper Functions
def get_user_input(prompt, default=None, is_boolean=False, choices=None):
    while True:
        user_input = (
            input(f"{prompt} [{'/'.join(choices) if choices else default}]: ")
            or default
        )
        if is_boolean:
            return user_input.lower() == "yes"
        if choices and user_input not in choices:
            print(f"Invalid choice. Please choose from {', '.join(choices)}.")
            continue
        return user_input


def setup_configuration(assay, genome, template_data):
    username = os.getenv("USER", "unknown_user")
    today = datetime.datetime.now().strftime("%Y-%m-%d")
    project_name = get_user_input(
        "What is your project name?", default=f"{username}_project"
    )
    project_name = project_name.replace(" ", "_")

    common_config = {
        "username": username,
        "project_date": today,
        "project_name": project_name,
        "genome": genome,
    }

    template_data.update(common_config)

    with open(os.path.join(template_dir, "preset_genomes.json"), "r") as f:
        genome_values = json.load(f)

    genome_dict = {}

    if genome == "other":
        genome = get_user_input("What is your genome name?", default="other")
        genome_dict = {
            genome: {
                "indices": (
                    get_user_input("Path to Bowtie2 genome indices:")
                    if assay in ["chip", "atac"]
                    else get_user_input("Path to STAR v2.7.10b genome indices:")
                ),
                "chromosome_sizes": get_user_input("Path to chromosome sizes file:"),
                "gtf": get_user_input("Path to GTF file:"),
                "blacklist": get_user_input("Path to blacklist bed file:"),
            }
        }
    else:
        if genome in genome_values:
            genome_dict[genome] = {
                "indices": genome_values[genome].get(
                    "bt2_indices" if assay in ["chip", "atac"] else "star_indices", ""
                ),
                "chromosome_sizes": genome_values[genome].get("chromosome_sizes", ""),
                "gtf": genome_values[genome].get("gtf", ""),
                "blacklist": genome_values[genome].get("blacklist", ""),
            }

    genome_config = {
        "genome": genome,
        "indices": genome_dict[genome]["indices"],
        "chromosome_sizes": genome_dict[genome]["chromosome_sizes"],
        "gtf": genome_dict[genome]["gtf"],
    }
    template_data.update(genome_config)

    # Fastqscreen
    template_data["fastq_screen"] = get_user_input(
        "Perform fastqscreen? (yes/no)", default="no", is_boolean=True
    )
    if template_data["fastq_screen"]:
        template_data["fastq_screen_config"] = get_user_input(
            "Path to fastqscreen config:",
            default="/ceph/project/milne_group/shared/seqnado_reference/fastqscreen_reference/fastq_screen.conf",
        )
    
    # Blacklist
    template_data["remove_blacklist"] = get_user_input(
        "Do you want to remove blacklist regions? (yes/no)",
        default="yes",
        is_boolean=True,
    )
    if template_data["remove_blacklist"]:
        template_data["blacklist"] = genome_dict[genome]["blacklist"]

    # Handle duplicates
    template_data["remove_pcr_duplicates"] = get_user_input(
        "Remove PCR duplicates? (yes/no)",
        default="yes" if assay in ["chip", "atac"] else "no",
        is_boolean=True,
    )
    if template_data["remove_pcr_duplicates"]:
        template_data["remove_pcr_duplicates_method"] = get_user_input(
            "Remove PCR duplicates method:", default="picard", choices=["picard"]
        )
        # Library Complexity
        template_data["library_complexity"] = get_user_input(
        "Calculate library complexity? (yes/no)", default="no", is_boolean=True
    )
    else:
        template_data["remove_pcr_duplicates_method"] = "False"
        template_data["library_complexity"] = "False"

    # Shift reads
    if assay == "atac":
        template_data["shift_atac_reads"] = (
            get_user_input(
                "Shift ATAC-seq reads? (yes/no)", default="yes", is_boolean=True
            )
            if assay == "atac"
            else "False"
        )

    # Spike in
    if assay in ["chip", "rna"]:
        template_data["spikein"] = get_user_input(
            "Do you have spikein? (yes/no)", default="no", is_boolean=True
        )
        if template_data["spikein"] and not assay == "rna":
            template_data["normalisation_method"] = get_user_input(
                "Normalisation method:",
                default="orlando",
                choices=["orlando", "with_input"],
            )
            template_data["reference_genome"] = get_user_input(
                "Reference genome:", default="hg38"
            )
            template_data["spikein_genome"] = get_user_input(
                "Spikein genome:", default="dm6"
            )

    # Make bigwigs
    template_data["make_bigwigs"] = get_user_input(
        "Do you want to make bigwigs? (yes/no)", default="no", is_boolean=True
    )
    if template_data["make_bigwigs"]:
        template_data["pileup_method"] = get_user_input(
            "Pileup method:", default="deeptools", choices=["deeptools", "homer"]
        )
        template_data["scale"] = get_user_input(
            "Scale bigwigs? (yes/no)", default="no", is_boolean=True
        )
        template_data["make_heatmaps"] = get_user_input(
            "Do you want to make heatmaps? (yes/no)", default="no", is_boolean=True
        )

    # Call peaks
    if assay in ["chip", "atac"]:
        template_data["call_peaks"] = get_user_input(
            "Do you want to call peaks? (yes/no)", default="no", is_boolean=True
        )
        if template_data["call_peaks"]:
            template_data["peak_calling_method"] = get_user_input(
                "Peak caller:",
                default="lanceotron",
                choices=["lanceotron", "macs", "homer"],
            )

    # Run DESeq2
    template_data["run_deseq2"] = (
        get_user_input("Run DESeq2? (yes/no)", default="no", is_boolean=True)
        if assay == "rna"
        else "False"
    )

    # Make UCSC hub
    template_data["make_ucsc_hub"] = get_user_input(
        "Do you want to make a UCSC hub? (yes/no)", default="no", is_boolean=True
    )

    template_data["UCSC_hub_directory"] = (
        get_user_input("UCSC hub directory:", default="/path/to/ucsc_hub/")
        if template_data["make_ucsc_hub"]
        else "."
    )
    template_data["email"] = (
        get_user_input("What is your email address?", default=f"{username}@example.com")
        if template_data["make_ucsc_hub"]
        else f"{username}@example.com"
    )
    template_data["color_by"] = (
        get_user_input("Color by (for UCSC hub):", default="samplename")
        if template_data["make_ucsc_hub"]
        else "samplename"
    )

    template_data["options"] = (
        TOOL_OPTIONS if assay in ["chip", "atac"] else TOOL_OPTIONS_RNA
    )


TOOL_OPTIONS = """
trim_galore:
    threads: 4
    options: --2colour 20 

bowtie2:
    threads: 8
    options:

picard:
    threads: 8
    options:

homer:
    use_input: true
    maketagdirectory:
    makebigwig:
    findpeaks:

deeptools:
    threads: 16
    alignmentsieve: --minMappingQuality 30 
    bamcoverage: --extendReads -bs 1 --normalizeUsing RPKM

macs:
    version: 2
    callpeak:

lanceotron:
    use_input: True
    callpeak: -c 0.5

heatmap:
    options:
    colormap: RdYlBu_r
"""

TOOL_OPTIONS_RNA = """
trim_galore:
    threads: 4
    options: --2colour 20 

star:
    threads: 16
    options: --quantMode TranscriptomeSAM GeneCounts --outSAMunmapped Within --outSAMattributes Standard

picard:
    threads: 8
    options:

featurecounts:
    threads: 16
    options: -s 0 -p --countReadPairs -t exon -g gene_id

homer:
    maketagdirectory:
    makebigwig:

deeptools:
    threads: 16
    alignmentsieve: --minMappingQuality 30 
    bamcoverage: -bs 1 --normalizeUsing CPM

heatmap:
    options:
    colormap: RdYlBu_r
"""


def create_config(assay, genome, rerun, debug=False):
    env = Environment(loader=FileSystemLoader(template_dir), auto_reload=False)

    template = env.get_template("config.yaml.jinja")
    template_deseq2 = env.get_template("deseq2.qmd.jinja")

    # Initialize template data
    template_data = {"assay": assay, "genome": genome}

    # Setup configuration
    setup_configuration(assay, genome, template_data)

    # Create directory and render template
    if rerun:
        dir_name = os.getcwd()
        with open(os.path.join(dir_name, f"config_{assay}.yml"), "w") as file:
            file.write(template.render(template_data))
    else:
        dir_name = f"{template_data['project_date']}_{template_data['assay']}_{template_data['project_name']}"
        os.makedirs(dir_name, exist_ok=True)
        fastq_dir = os.path.join(dir_name, "fastq")
        os.makedirs(fastq_dir, exist_ok=True)

        with open(os.path.join(dir_name, f"config_{assay}.yml"), "w") as file:
            file.write(template.render(template_data))

    # add deseq2 qmd file if rna
    if assay == "rna":
        with open(
            os.path.join(dir_name, f"deseq2_{template_data['project_name']}.qmd"), "w"
        ) as file:
            file.write(template_deseq2.render(template_data))

    print(
        f"Directory '{dir_name}' has been created with the 'config_{assay}.yml' file."
    )

    if debug:
        with open(os.path.join(dir_name, "data.json"), "w") as file:
            json.dump(template_data, file)
