import os
import csv
from pathlib import Path
import pandas as pd
from seqnado.inputs import Assay, FastqCollection, FastqCollectionForIP
from seqnado.inputs.validation import DesignDataFrame

def get_designs_dir():
    """Returns the directory where designs are stored (default local designs/ folder)."""
    designs_dir = Path("designs")
    designs_dir.mkdir(exist_ok=True)
    return designs_dir

def list_designs(directory: str = None):
    """
    List design files in the given directory (or default designs dir).
    Returns a list of dictionaries with 'name', 'path', 'modified'.
    """
    if directory:
        target_dir = Path(directory)
    else:
        target_dir = get_designs_dir()
        
    if not target_dir.exists():
        return []

    designs = []
    # Only list .csv files in the top level if it's the default dir
    # If it's a scan, we might want recursive? Let's check user intent.
    # User said "identifying them in folders". 
    # Let's simple listing for the default view.
    
    try:
        for f in target_dir.glob("*.csv"):
            designs.append({
                "name": f.name,
                "path": str(f.resolve()),
                "modified": os.path.getmtime(f)
            })
    except Exception:
        pass # Handle permission errors etc
        
    return designs

def scan_for_designs(root_path: str, max_depth: int = 3):
    """
    Scans a directory recursively for CSV files that might be designs.
    """
    root = Path(root_path).resolve()
    found_designs = []
    
    # Simple walk with depth limit
    initial_depth = len(root.parts)
    
    for dirpath, dirnames, filenames in os.walk(root):
        current_depth = len(Path(dirpath).parts)
        if current_depth - initial_depth > max_depth:
            del dirnames[:] # Don't go deeper
            continue
            
        for f in filenames:
            if f.endswith(".csv"):
                full_path = Path(dirpath) / f
                found_designs.append({
                    "name": f,
                    "rel_path": str(full_path.relative_to(root)),
                    "path": str(full_path),
                    "modified": os.path.getmtime(full_path)
                })
                
    return found_designs

def load_design(filename_or_path: str):
    """
    Loads a design CSV. Accepts a simple filename (assumed in designs/) OR an absolute path.
    """
    # Check if it's an absolute path or exists locally
    p = Path(filename_or_path)
    if p.is_absolute():
        path = p
    elif p.exists():
        path = p
    # Check if it was an absolute path but got stripped (e.g. by Flask routing)
    elif Path(f"/{filename_or_path}").exists():
        path = Path(f"/{filename_or_path}")
    else:
        path = get_designs_dir() / filename_or_path

    if not path.exists():
        # Last ditch effort: maybe it IS absolute but doesn't exist yet (unlikely for load, but logic holds)
        if Path(f"/{filename_or_path}").parent.exists():
             path = Path(f"/{filename_or_path}")
        else:
             return []
    
    try:
        # data = pd.read_csv(path).to_dict(orient="records") # Pandas is robust but we want to preserve structure
        # Just simple dict reader for now to pass to JS
        with open(path, "r") as f:
            reader = csv.DictReader(f)
            data = list(reader)
        # Store strict filename/path for save logic?
        # Actually save_design receives the same filename argument from the frontend usually.
        return data
    except Exception as e:
        print(f"Error loading design: {e}")
        return []

def save_design(filename_or_path: str, data: list):
    """
    Saves a design CSV. Accepts name (saves to designs/) or absolute path.
    """
    p = Path(filename_or_path)
    # If it looks like a path (has separators) or is absolute, use it.
    if p.is_absolute():
         path = p
    elif Path(f"/{filename_or_path}").parent.exists(): # Heuristic for stripped absolute path
         path = Path(f"/{filename_or_path}")
    elif len(p.parts) > 1:
         # Relative path provided?
         path = p
    else:
         path = get_designs_dir() / filename_or_path
         
    if not data:
        # Create empty file
        try:
            path.touch()
            return True
        except Exception as e:
            print(f"Error saving design: {e}")
            return False

    keys = data[0].keys()
    
    try:
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(data)
        return True
    except Exception as e:
        print(f"Error saving design: {e}")
        return False

def create_new_design(filename: str):
    """Creates a new empty design in the default directory."""
    if not filename.endswith(".csv"):
        filename += ".csv"
    path = get_designs_dir() / filename
    # Create with basic columns? Or empty. 
    # Let's create with 'sample_id' as it's required for everything
    with open(path, "w", newline="") as f:
        f.write("sample_id,fastq_1,fastq_2,condition\n")
    return str(path)

def apply_schema_defaults(df: pd.DataFrame, assay: Assay) -> pd.DataFrame:
    """
    Applies schema defaults to the DataFrame, adding missing columns.
    Mimics CLI behavior but non-interactively.
    """
    try:
        schema = DesignDataFrame.to_schema()
    except Exception:
        return df

    to_ignore = {"r2", "r1_control", "r2_control", "ip", "control"}
    
    # Iterate over schema columns
    for name, col in schema.columns.items():
        if name in df.columns:
            continue
        if name in to_ignore:
            continue
            
        default = getattr(col, "default", None)
        nullable = bool(getattr(col, "nullable", False))
        
        if default is not None:
            df[name] = default
        elif nullable:
            df[name] = "" # Fill with empty string for easier editing
            
    return df

def generate_design_from_fastqs(folder_path: str, assay_name: str, output_name: str) -> str:
    """
    Generates a design CSV from FASTQ files in a given folder.
    
    Args:
        folder_path: Directory containing .fastq.gz files.
        assay_name: Name of the assay (e.g., 'rna', 'atac').
        output_name: Filename for the output CSV.
    
    Returns:
        Path to the saved CSV file.
    """
    folder = Path(folder_path).resolve()
    if not folder.exists():
         raise FileNotFoundError(f"Folder not found: {folder_path}")

    try:
        assay = Assay.from_clean_name(assay_name)
    except ValueError:
        raise ValueError(f"Invalid assay: {assay_name}")

    # Use SeqNado input logic with built-in discovery (heuristics)
    if assay in {Assay.CHIP, Assay.CAT}:
        design_obj = FastqCollectionForIP.from_directory(
            assay=assay,
            directory=folder
        )
    else:
        design_obj = FastqCollection.from_directory(
            assay=assay,
            directory=folder
        )
        
    df = design_obj.to_dataframe().sort_values("sample_id")
    
    # Apply heuristic schema defaults
    df = apply_schema_defaults(df, assay)
    
    # Save to designs dir (or specified output path if absolute??)
    # User usually wants it in 'designs/'
    if not output_name.endswith(".csv"):
        output_name += ".csv"
        
    # If output_name is just a name, go to designs_dir. If absolute, use it.
    out_p = Path(output_name)
    if not out_p.is_absolute():
        final_path = get_designs_dir() / output_name
    else:
        final_path = out_p
        
    final_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(final_path, index=False)
    
    return str(final_path)
