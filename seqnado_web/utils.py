import os
import json
from pathlib import Path
from typing import Dict, Any, List

from seqnado import Assay
from seqnado.config import load_genome_configs, GenomeConfig

def get_genome_config_path() -> Path:
    """Get the path to the genome configuration file."""
    return Path(os.getenv("SEQNADO_CONFIG", Path.home())) / ".config/seqnado/genome_config.json"

def load_all_genomes() -> Dict[str, Any]:
    """
    Load all genome configurations. 
    Returns a dictionary of raw config data (dicts), not GenomeConfig objects,
    to be easily serializable for templates.
    """
    config_path = get_genome_config_path()
    if not config_path.exists():
        return {}
    
    with open(config_path, "r") as f:
        return json.load(f)

def save_genome_config(name: str, config_data: Dict[str, Any]) -> None:
    """Save or update a genome configuration."""
    all_configs = load_all_genomes()
    all_configs[name] = config_data
    
    config_path = get_genome_config_path()
    config_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(config_path, "w") as f:
        json.dump(all_configs, f, indent=4)

def delete_genome_config(name: str) -> None:
    """Delete a genome configuration."""
    all_configs = load_all_genomes()
    if name in all_configs:
        del all_configs[name]
        
        config_path = get_genome_config_path()
        with open(config_path, "w") as f:
            json.dump(all_configs, f, indent=4)

def get_available_assays() -> List[str]:
    """Get list of available assays."""
    return Assay.all_assay_clean_names()
