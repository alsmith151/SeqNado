import os
from pathlib import Path
import socket
from flask import Flask, render_template, request, redirect, url_for, flash
from seqnado_web.utils import load_all_genomes, save_genome_config, delete_genome_config

app = Flask(__name__)
app.secret_key = os.urandom(24)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/genomes")
def genomes():
    genomes_data = load_all_genomes()
    return render_template("genomes.html", genomes=genomes_data)

@app.route("/genomes/add", methods=["GET", "POST"])
def add_genome():
    if request.method == "POST":
        name = request.form.get("name")
        if not name:
            # Handle error
            return redirect(url_for("add_genome"))
            
        organism = request.form.get("organism")
        if not organism:
            # Simple heuristic prediction mimicking GenomeConfig.predict_organism
            if "hg" in name:
                organism = "Homo sapiens"
            elif "mm" in name:
                organism = "Mus musculus"
            else:
                organism = "Unknown"

        config_data = {
            "name": name,
            "organism": organism,
            "version": request.form.get("version"),
            "fasta": request.form.get("fasta") or None,
            "gtf": request.form.get("gtf") or None,
            "star_index": request.form.get("star_index") or None,
            "bt2_index": request.form.get("bt2_index") or None,
            "blacklist": request.form.get("blacklist") or None,
        }
        # Clean None values
        config_data = {k: v for k, v in config_data.items() if v is not None}
        
        save_genome_config(name, config_data)
        return redirect(url_for("genomes"))
        
    return render_template("genome_form.html", is_edit=False, name="", config={})

# ... (omitted handlers) ...

@app.route("/api/browse", methods=["GET"])
def api_browse():
    path = request.args.get("path", ".")
    
    # Security/Safety check: maybe restrict to certain root? 
    # For a dev tool running on HPC/local as user, usually we allow full access like a terminal.
    # We resolve it to absolute path.
    if path == ".":
        abs_path = Path.cwd()
    else:
        abs_path = Path(path).resolve()
        
    if not abs_path.exists():
         return {"error": "Path does not exist"}, 400
         
    # If it's a file, return parent? Or handle it?
    # Usually browsing implies listing a directory. 
    if abs_path.is_file():
         abs_path = abs_path.parent
         
    try:
        items = []
        # Parent directory option (unless root)
        if abs_path.parent != abs_path:
             items.append({
                "name": "..",
                "path": str(abs_path.parent),
                "is_dir": True
            })
            
        for p in abs_path.iterdir():
            try:
                # Show both dirs and files, but maybe filter out hidden?
                if not p.name.startswith("."): 
                    items.append({
                        "name": p.name,
                        "path": str(p),
                        "is_dir": p.is_dir()
                    })
            except PermissionError:
                continue
                
        # Sort items: directories first, then files
        items.sort(key=lambda x: (not x["is_dir"], x["name"]))
        
        return {
            "current_path": str(abs_path),
            "items": items
        }
    except Exception as e:
        return {"error": str(e)}, 500

@app.route("/genomes/<name>", methods=["GET", "POST"])
def edit_genome(name):
    if request.method == "POST":
        action = request.form.get("action")
        if action == "delete":
            delete_genome_config(name)
            return redirect(url_for("genomes"))
            
        organism = request.form.get("organism")
        if not organism:
            # Simple heuristic prediction mimicking GenomeConfig.predict_organism
            if "hg" in name:
                organism = "Homo sapiens"
            elif "mm" in name:
                organism = "Mus musculus"
            else:
                organism = "Unknown"

        config_data = {
            "name": name,
            "organism": organism,
            "version": request.form.get("version"),
            "fasta": request.form.get("fasta") or None,
            "gtf": request.form.get("gtf") or None,
            "star_index": request.form.get("star_index") or None,
            "bt2_index": request.form.get("bt2_index") or None,
            "blacklist": request.form.get("blacklist") or None,
        }
        # Clean None values
        config_data = {k: v for k, v in config_data.items() if v is not None}
        
        save_genome_config(name, config_data)
        return redirect(url_for("genomes"))

    genomes_data = load_all_genomes()
    config = genomes_data.get(name)
    if not config:
        return redirect(url_for("genomes"))
        
    return render_template("genome_form.html", is_edit=True, name=name, config=config)

@app.route("/configs")
def configs():
    return render_template("index.html") # Placeholder

@app.route("/run")
def run_experiment():
    return render_template("run.html") # Simplified

@app.route("/init")
def init_project():
    from seqnado_web.utils import get_available_assays, load_all_genomes
    assays = get_available_assays()
    genomes = load_all_genomes()
    return render_template("init.html", assays=assays, genomes=genomes)

@app.route("/init/scan-fastqs", methods=["POST"])
def scan_fastqs_api():
    from seqnado_web.init_utils import scan_directory_for_fastqs
    path = request.form.get("path")
    if not path:
        return {"error": "Path required"}, 400
    files = scan_directory_for_fastqs(path)
    return {"files": files}

@app.route("/init/create", methods=["POST"])
def create_project():
    from seqnado_web.init_utils import create_project_structure
    
    assay = request.form.get("assay")
    project_name = request.form.get("project_name")
    genome = request.form.get("genome")
    output_dir = request.form.get("output_dir", ".")
    
    selected_fastqs = request.form.getlist("selected_fastqs")
    
    # Gather config options from Step 4 form
    config_options = {
        "make_bigwigs": request.form.get("make_bigwigs") == "yes",
        "pileup_method": request.form.get("pileup_method", "deeptools"),
        "binsize": int(request.form.get("binsize", 10)),
        "make_heatmaps": request.form.get("make_heatmaps") == "yes",
        "make_ucsc_hub": request.form.get("make_ucsc_hub") == "yes",
        # Assay-specific
        "tn5_shift": request.form.get("tn5_shift") == "yes",
        "call_peaks": request.form.get("call_peaks") == "yes",
        "peak_calling_method": request.form.get("peak_calling_method", "lanceotron"),
        "has_spikein": request.form.get("has_spikein") == "yes",
        "quant_method": request.form.get("quant_method", "feature_counts"),
        "run_deseq2": request.form.get("run_deseq2") == "yes",
        "mcc_viewpoints": request.form.get("mcc_viewpoints"),
        "mcc_resolutions": request.form.get("mcc_resolutions", "100,1000"),
    }
    
    print(f"DEBUG: create_project called with:")
    print(f"  Assay: {assay}")
    print(f"  Project: {project_name}")
    print(f"  Output: {output_dir}")
    print(f"  Selected Fastqs ({len(selected_fastqs)}): {selected_fastqs[:3]}...")
    print(f"  Config Options: {config_options}")
    
    if not all([assay, project_name, genome]):
        flash("Missing required fields", "error")
        return redirect(url_for("init_project"))
        
    try:
        run_dir = create_project_structure(
            assay_name=assay,
            project_name=project_name,
            genome_name=genome,
            output_dir=output_dir,
            selected_fastqs=selected_fastqs,
            config_options=config_options
        )
        flash(f"Project created successfully at {run_dir}", "success")
        # Go to designs list? Or directly to the project runner?
        return redirect(url_for("run_experiment"))
    except Exception as e:
        import traceback
        traceback.print_exc()
        flash(f"Error creating project: {str(e)}", "error")
        return redirect(url_for("init_project"))

@app.route("/api/load_project_config", methods=["POST"])
def api_load_project_config():
    import yaml
    path = request.form.get("path")
    if not path or not os.path.exists(path):
        return {"error": "Invalid path"}, 400
    
    # Check for config_*.yaml
    try:
        path_p = Path(path)
        configs = list(path_p.glob("config_*.yaml"))
        if not configs:
            return {"error": "No config_*.yaml found in directory"}, 404
            
        config_path = configs[0]
        # Parse basic details
        with open(config_path) as f:
            cfg = yaml.safe_load(f)
            
        return {
            "assay": cfg.get("assay") if isinstance(cfg.get("assay"), str) else str(cfg.get("assay", "?")), # enum handling?
            "project": cfg.get("project", {}).get("name", "Unknown"),
            "genome": cfg.get("genome", {}).get("name", "Unknown"),
            "config_file": config_path.name
        }
    except Exception as e:
        return {"error": str(e)}, 500

@app.route("/run/submit", methods=["POST"])
def submit_run():
    project_dir = request.form.get("project_dir")
    profile = request.form.get("profile")
    unlock = request.form.get("unlock") == "yes"
    
    if not project_dir:
        return "Project directory required", 400

    try:
        from seqnado_web.run_utils import get_run_command
        
        # Find config again
        configs = list(Path(project_dir).glob("config_*.yaml"))
        if not configs:
            return "Config not found", 404
        config_path = configs[0]
        
        run_command = get_run_command(config_path, profile=profile)
        
        if unlock:
            # Prepend unlock command
            unlock_cmd = run_command + " --unlock"
            # run_command = f"{unlock_cmd} && {run_command}" # Display both?
            # Or just pass the unlock flag logic separately.
            # Let's just output instructions for now or execute.
            pass

        # Since we are verified execution, let's just show it for safety as requested in previous turn
        # "Run Snakemake pipeline immediately..." - The form said "Run Pipeline".
        # Let's assume we WANT to run it or at least show the command to run.
        
        run_now = True # The form action implies running or preparing to run.
        
        return render_template("run_result.html", 
                               project_name=Path(project_dir).name, 
                               run_command=run_command,
                               run_dir=project_dir,
                               run_now=run_now,
                               unlock=unlock)
        
    except Exception as e:
         return f"Error: {str(e)}", 500

# --- Design Editor Routes ---
from seqnado_web.design_utils import list_designs, load_design, save_design, create_new_design, get_designs_dir

@app.route("/designs")
def designs_list():
    from seqnado_web.utils import get_available_assays
    
    # Ensure designs dir exists
    get_designs_dir()
    designs = list_designs()
    assays = get_available_assays()
    return render_template("designs.html", designs=designs, assays=assays, scan_mode=False)

@app.route("/designs/create", methods=["POST"])
def create_design():
    filename = request.form.get("filename")
    if filename:
        create_new_design(filename)
    return redirect(url_for("designs_list"))

@app.route("/designs/<path:filename>")
def edit_design(filename):
    data = load_design(filename)
    if data is None:
        return redirect(url_for("designs_list"))
    return render_template("design_editor.html", filename=filename, data=data)

@app.route("/designs/scan", methods=["GET", "POST"])
def scan_designs():
    from seqnado_web.design_utils import scan_for_designs
    from seqnado_web.utils import get_available_assays
    
    if request.method == "POST":
        path = request.form.get("scan_path")
        if not path or not os.path.exists(path):
            flash("Path does not exist", "error")
            return redirect(url_for("designs_list")) 
            
        found_designs = scan_for_designs(path)
        assays = get_available_assays()
        return render_template("designs.html", designs=found_designs, assays=assays, scan_path=path, scan_mode=True)
    
    return redirect(url_for("designs_list"))

@app.route("/designs/generate", methods=["POST"])
def generate_design():
    from seqnado_web.design_utils import generate_design_from_fastqs
    
    folder_path = request.form.get("folder_path")
    assay_name = request.form.get("assay")
    output_name = request.form.get("output_name")
    output_dir = request.form.get("output_dir")
    
    if not all([folder_path, assay_name, output_name]):
        flash("All fields required", "error")
        return redirect(url_for("designs_list"))
        
    try:
        # Handle custom output directory
        if output_dir:
            if not output_name.endswith(".csv"):
                 output_name += ".csv"
            # Combine dir and name into an absolute path or relative path treated as path
            full_output_path = str(Path(output_dir) / output_name)
        else:
            full_output_path = output_name # Let utils handle default logic

        path = generate_design_from_fastqs(folder_path, assay_name, full_output_path)
        flash(f"Design generated successfully: {path}", "success")
        return redirect(url_for("edit_design", filename=path))
    except Exception as e:
        flash(f"Error generating design: {str(e)}", "error")
        return redirect(url_for("designs_list"))

@app.route("/api/design/save", methods=["POST"], endpoint="api_save_design")
def api_save_design():
    filename = request.json.get("filename")
    data = request.json.get("data")
    if save_design(filename, data):
        return {"status": "success"}, 200
    else:
        return {"status": "error"}, 500


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5001)
