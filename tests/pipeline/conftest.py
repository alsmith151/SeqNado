"""Heavyweight pipeline fixtures for Snakemake E2E tests.

These fixtures download small reference bundles into tests/data on first use,
prepare a temporary run directory, and configure environment variables so that
the SeqNado CLI and Snakemake pipeline can run under pytest.

They are intentionally scoped at session/function levels to balance reuse and
isolation across assays when --assays is used.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
from datetime import datetime

import pytest
import requests


# ------------------
# Path convenience
# ------------------

@pytest.fixture(scope="session")
def repo_path() -> Path:
    return Path(__file__).resolve().parents[2]


@pytest.fixture(scope="session")
def package_path(repo_path: Path) -> Path:
    return repo_path / "seqnado"


@pytest.fixture(scope="session")
def test_dir_path(repo_path: Path) -> Path:
    return repo_path / "tests"


@pytest.fixture(scope="session")
def test_data_path(test_dir_path: Path) -> Path:
    p = test_dir_path / "data"
    p.mkdir(parents=True, exist_ok=True)
    return p


@pytest.fixture(scope="session")
def workflow_path(package_path: Path) -> Path:
    return package_path / "workflow"


@pytest.fixture(scope="session")
def test_profile_path(workflow_path: Path) -> Path:
    return workflow_path / "envs" / "profiles" / "profile_test"


@pytest.fixture(scope="session")
def genome_path(test_data_path: Path) -> Path:
    p = test_data_path / "genome"
    p.mkdir(parents=True, exist_ok=True)
    return p


# ------------------
# Reference bundles
# ------------------

@pytest.fixture(scope="function")
def genome_index_path(genome_path: Path, assay: str) -> Path:
    if "rna" in assay:
        return genome_path / "STAR_chr21_rna_spikein"
    elif "meth" in assay:
        return genome_path / "bt2_chr21_meth"
    else:
        return genome_path / "bt2_chr21_dm6_chr2L"


@pytest.fixture(scope="function")
def index(genome_index_path: Path, genome_path: Path) -> Path:
    download_index = not genome_index_path.exists()
    suffix = genome_index_path.with_suffix(".tar.gz").name
    url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/{suffix}"

    if "STAR" in str(genome_index_path):
        indicies_path = genome_index_path
    elif "meth" in str(genome_index_path):
        indicies_path = genome_index_path / "chr21_meth"
    else:
        indicies_path = genome_index_path / "bt2_chr21_dm6_chr2L"

    if download_index:
        r = requests.get(url, stream=True)
        r.raise_for_status()
        tar_index = genome_index_path.with_suffix(".tar.gz")
        with open(tar_index, "wb") as f:
            f.write(r.content)
        with tarfile.open(tar_index) as tar:
            if "bt2" in str(genome_index_path):
                genome_index_path.mkdir(parents=True, exist_ok=True)
                tar.extractall(path=genome_index_path, filter="data")
            else:
                tar.extractall(path=genome_path, filter="data")
        os.remove(tar_index)

    return indicies_path


@pytest.fixture(scope="session")
def chromsizes(genome_path: Path) -> Path:
    suffix = "chr21.fa.fai"
    dest = genome_path / suffix

    if not (genome_path / "chr21_rename.fa.fai").exists():
        url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/{suffix}"
        r = requests.get(url, stream=True)
        r.raise_for_status()
        with open(dest, "wb") as f:
            f.write(r.content)

            # Append dm6 chromsizes
            url2 = "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes"
            r2 = requests.get(url2, stream=True)
            r2.raise_for_status()
            for line in r2.iter_lines():
                f.write(b"dm6_" + line + b"\n")

    return dest


@pytest.fixture(scope="function")
def gtf(genome_path: Path, assay: str, index: Path) -> Path:
    if "rna" in assay:
        gtf_path = genome_path / "chr21_rna_spikein.gtf"
    else:
        gtf_path = genome_path / "chr21.gtf"

    if not gtf_path.exists():
        url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/{gtf_path.name}"
        r = requests.get(url, stream=True)
        r.raise_for_status()
        with open(gtf_path, "wb") as f:
            f.write(r.content)
    return gtf_path


@pytest.fixture(scope="function")
def blacklist(genome_path: Path) -> Path:
    blacklist_path = genome_path / "hg38-blacklist.v2.bed.gz"
    if not blacklist_path.exists():
        url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
        r = requests.get(url, stream=True)
        r.raise_for_status()
        with open(blacklist_path, "wb") as f:
            f.write(r.content)
    return blacklist_path


@pytest.fixture(scope="session")
def genome_files(genome_path: Path) -> tuple[Path, Path]:
    fasta = genome_path / "chr21_meth.fa"
    fasta_fai = genome_path / "chr21_meth.fa.fai"
    if not fasta.exists():
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/chr21_meth.fa"
        r = requests.get(url, stream=True)
        r.raise_for_status()
        with open(fasta, "wb") as f:
            f.write(r.content)
    if not fasta_fai.exists():
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/chr21_meth.fa.fai"
        r = requests.get(url, stream=True)
        r.raise_for_status()
        with open(fasta_fai, "wb") as f:
            f.write(r.content)
    return fasta, fasta_fai


@pytest.fixture(scope="session")
def mcc_files(genome_path: Path) -> dict[str, Path]:
    files: dict[str, Path] = {}
    url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/seqnado_data/test_viewpoints.bed"
    r = requests.get(url, stream=True)
    r.raise_for_status()
    mcc_viewpoints = genome_path / "mcc_viewpoints.bed"
    with open(mcc_viewpoints, "wb") as f:
        f.write(r.content)
    files["viewpoints"] = mcc_viewpoints
    return files


# ------------------
# Assay helpers
# ------------------

@pytest.fixture(scope="function")
def assay_type(assay: str) -> str:
    return re.sub(r"(.*)\-.*", r"\1", assay)


@pytest.fixture(scope="session")
def fastqs(test_data_path: Path, selected_assays: list[str]) -> dict[str, list[Path]]:
    """Download once and index by assay; each test copies needed files."""
    target_dir = test_data_path / "fastq"
    target_dir.mkdir(parents=True, exist_ok=True)

    # Download archive only if directory is empty
    if not any(target_dir.glob("*.fastq.gz")):
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/seqnado_data/fastq.tar.gz"
        r = requests.get(url, stream=True)
        r.raise_for_status()
        tar_path = target_dir.parent / "fastq.tar.gz"
        with open(tar_path, "wb") as f:
            f.write(r.content)
        temp_extract_dir = target_dir.parent / "temp_fastq_extracted"
        temp_extract_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(path=temp_extract_dir)
        for fastq_file in temp_extract_dir.rglob("*.fastq.gz"):
            shutil.move(str(fastq_file), target_dir / fastq_file.name)
        shutil.rmtree(temp_extract_dir)
        tar_path.unlink(missing_ok=True)

    def pick(assay_name: str) -> list[Path]:
        m = {
            "atac": "atac*.fastq.gz",
            "chip": "chip-rx_*.fastq.gz",
            "chip-rx": "chip-rx_*.fastq.gz",
            "rna": "rna_*.fastq.gz",
            "rna-rx": "rna-spikein*.fastq.gz",
            "snp": "snp*.fastq.gz",
            "cat": "chip-rx_*.fastq.gz",
            "meth": "meth*.fastq.gz",
            "mcc": "mcc*.fastq.gz",
        }
        pattern = m.get(assay_name)
        if not pattern:
            raise ValueError(f"Unsupported assay: {assay_name}")
        return sorted(target_dir.glob(pattern))

    return {a: pick(a) for a in selected_assays}


@pytest.fixture(scope="function")
def run_directory(tmp_path_factory: pytest.TempPathFactory, assay: str) -> Path:
    return tmp_path_factory.mktemp(assay)


@pytest.fixture(scope="function")
def plot_bed(test_data_path: Path) -> Path:
    return test_data_path / "plotting_coordinates.bed"


@pytest.fixture(scope="function", autouse=True)
def set_up_run_dir(run_directory: Path, fastqs: dict[str, list[Path]], assay: str):
    """Chdir to run dir and copy fastqs needed for this assay."""
    cwd = Path.cwd()
    os.chdir(run_directory)
    try:
        for fq in fastqs[assay]:
            shutil.copy(fq, run_directory)
        yield
    finally:
        os.chdir(cwd)


@pytest.fixture(scope="function", autouse=True)
def run_init(index: Path, chromsizes: Path, gtf: Path, blacklist: Path, run_directory: Path, assay: str, monkeypatch: pytest.MonkeyPatch, genome_files: tuple[Path, Path], genome_path: Path):
    """Run `seqnado init` and materialize genome_config.json with appropriate indices."""
    monkeypatch.setenv("SEQNADO_CONFIG", str(run_directory))
    monkeypatch.setenv("HOME", str(run_directory))

    process = subprocess.run(["seqnado", "init"], input="y\n", text=True, cwd=run_directory, capture_output=True)

    genome_config_file = run_directory / ".config" / "seqnado" / "genome_config.json"
    if "rna" in assay:
        index_dict = {"star_index": str(index), "bt2_index": "NA"}
    else:
        index_dict = {"star_index": "NA", "bt2_index": str(index)}

    fasta, _ = genome_files
    genome_config_dict = {
        "hg38": {
            "star_index": index_dict["star_index"],
            "bt2_index": index_dict["bt2_index"],
            "chromosome_sizes": str(chromsizes),
            "gtf": str(gtf),
            "blacklist": str(blacklist),
            "genes": "",
            "fasta": str(fasta) if assay in ["meth", "snp"] else "",
        }
    }
    genome_config_file.parent.mkdir(parents=True, exist_ok=True)
    with open(genome_config_file, "w") as f:
        json.dump(genome_config_dict, f, indent=4)

    # Copy fastq_screen.conf from test data, updating the index path
    source_fastq_screen_config = genome_path / "fastq_screen.conf"
    dest_fastq_screen_config = run_directory / ".config" / "seqnado" / "fastq_screen.conf"
    dest_fastq_screen_config.parent.mkdir(parents=True, exist_ok=True)
    
    if source_fastq_screen_config.exists():
        # Read template and replace index path
        with open(source_fastq_screen_config, "r") as f:
            config_content = f.read()
        # Replace the template path with actual index path
        config_content = config_content.replace(
            "/path",
            str(index)
        )
        with open(dest_fastq_screen_config, "w") as f:
            f.write(config_content)
    else:
        # Fallback: create minimal config
        with open(dest_fastq_screen_config, "w") as f:
            f.write(f"""# Test fastq_screen.conf file

DATABASE\tTest\t{index}
""")

    assert process.returncode == 0, f"seqnado init failed: {process.stderr}"
    assert genome_config_file.exists(), "genome_config.json missing after init"


@pytest.fixture(scope="function")
def config_yaml(run_directory: Path, assay_type: str, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.setenv("SEQNADO_CONFIG", str(run_directory))
    monkeypatch.setenv("HOME", str(run_directory))

    # Create seqnado_output dir for UCSCHubConfig validation
    (run_directory / "seqnado_output").mkdir(exist_ok=True)

    # Use subprocess with --no-interactive and --render-options for non-interactive config generation
    date = datetime.now().strftime("%Y-%m-%d")
    project_dir = run_directory / f"{date}_{assay_type}_test"
    project_dir.mkdir(parents=True, exist_ok=True)
    config_file_path = project_dir / f"config_{assay_type}.yaml"
    
    result = subprocess.run(
        [
            "seqnado",
            "config",
            assay_type,
            "--no-interactive",
            "--no-make-dirs",
            "--render-options",
            "-o",
            str(config_file_path),
        ],
        cwd=run_directory,
        capture_output=True,
        text=True,
    )
    
    assert result.returncode == 0, f"seqnado config failed: {result.stderr}\n{result.stdout}"
    assert config_file_path.exists(), f"{assay_type} config file not created."
    return config_file_path


@pytest.fixture(scope="function")
def config_yaml_for_testing(config_yaml: Path, assay: str) -> Path:
    import yaml
    with open(config_yaml, "r") as f:
        config = yaml.safe_load(f)

    match assay:
        case "chip":
            config["scale"] = "yes"
            config["library_complexity"] = False
        case "chip-rx":
            config["peak_calling_method"] = ["seacr"]
        case "atac":
            config["pileup_method"] = ["deeptools", "homer"]
            config["call_peaks"] = True
            config["peak_calling_method"] = ["lanceotron", "macs", "homer"]
        case "cat":
            config["pileup_method"] = ["deeptools", "bamnado"]
        case "mcc":
            config["call_peaks"] = False
            config["peak_calling_method"] = ["lanceotron-mcc"]
        case _:
            pass

    with open(config_yaml, "w") as f:
        yaml.dump(config, f)
    return Path(config_yaml)


@pytest.fixture(scope="function")
def seqnado_run_dir(config_yaml_for_testing: Path) -> Path:
    return Path(config_yaml_for_testing).parent


@pytest.fixture(scope="function")
def design(seqnado_run_dir: Path, assay_type: str, assay: str, run_directory: Path) -> Path:
    import pandas as pd
    
    # Copy FASTQs from run_directory to seqnado_run_dir
    for fq in run_directory.glob("*.fastq.gz"):
        shutil.copy(fq, seqnado_run_dir)
    
    metadata_file = f"metadata_{assay_type}.csv"
    cmd = ["seqnado", "design", assay_type, "-o", metadata_file, "--no-interactive", "--accept-all-defaults"]
    completed = subprocess.run(cmd, cwd=seqnado_run_dir, capture_output=True, text=True)
    assert completed.returncode == 0, f"seqnado design failed: {completed.stderr}"

    if assay == "chip":
        df = pd.read_csv(seqnado_run_dir / metadata_file, index_col=0)
        df["merge"] = "MLL-MERGED-TOGETHER"
        df.to_csv(seqnado_run_dir / metadata_file)
    elif assay == "rna-rx":
        df = pd.read_csv(seqnado_run_dir / metadata_file, index_col=0)
        df["deseq2"] = df.index.str.split("-").str[-2]
        df.to_csv(seqnado_run_dir / metadata_file)
    return seqnado_run_dir / metadata_file


@pytest.fixture(scope="function", autouse=True)
def apptainer_args(index: Path, test_data_path: Path, monkeypatch: pytest.MonkeyPatch):
    import importlib.resources
    import seqnado.data

    indicies_mount = index.parent if not index.is_dir() else index
    tmpdir = Path(os.environ.get("TMPDIR", "/tmp") or "/tmp")
    wd = Path(os.getcwd()).resolve()
    apptainer_cache_dir = Path.home() / ".apptainer" / "cache"
    multiqc_config = Path(importlib.resources.files(seqnado.data) / "multiqc_config.yaml").absolute().resolve()
    multiqc_config_parent = multiqc_config.parent

    monkeypatch.setenv(
        "APPTAINER_BINDPATH",
        f"{wd}:{wd},{test_data_path}:{test_data_path},{indicies_mount}:{indicies_mount},{tmpdir}:{tmpdir},{multiqc_config_parent}:{multiqc_config_parent}",
    )
    monkeypatch.setenv("APPTAINER_CACHEDIR", str(apptainer_cache_dir), prepend=False)
    monkeypatch.setenv("APPTAINER_TMPDIR", str(tmpdir), prepend=False)
