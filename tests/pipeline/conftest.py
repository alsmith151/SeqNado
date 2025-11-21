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
import re
import shutil
import subprocess
import tarfile
import time
from datetime import datetime
from pathlib import Path

import pytest
import requests


def _parse_assays(config: pytest.Config) -> list[str]:
    """Parse --assays option into list. Single source of truth for assay parsing."""
    assays_str: str = config.getoption("--assays") or "chip"
    assays = [a.strip() for a in assays_str.split(",") if a.strip()]
    return assays or ["chip"]


def _download_with_retry(
    url: str, dest: Path, max_retries: int = 3, timeout: int = 30
) -> None:
    """Download a file with retry logic."""
    for attempt in range(max_retries):
        try:
            r = requests.get(url, stream=True, timeout=timeout)
            r.raise_for_status()
            with open(dest, "wb") as f:
                f.write(r.content)
            return
        except (requests.RequestException, OSError) as e:
            if attempt < max_retries - 1:
                time.sleep(2**attempt)  # Exponential backoff
                continue
            raise RuntimeError(
                f"Failed to download {url} after {max_retries} attempts"
            ) from e


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
def star_index(genome_path: Path) -> Path:
    suffix = "STAR_chr21_rna_spikein.tar.gz"
    dest = genome_path / "STAR_chr21_rna_spikein"

    if not dest.exists():
        url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/{suffix}"
        _download_with_retry(url, genome_path / suffix)
        with tarfile.open(genome_path / suffix) as tar:
            tar.extractall(path=genome_path)
        os.remove(genome_path / suffix)

    return dest

@pytest.fixture(scope="session")
def bt2_index(genome_path: Path) -> Path:
    suffix = "bt2_chr21_dm6_chr2L.tar.gz"
    dest = genome_path / "bt2_chr21_dm6_chr2L/bt2_chr21_dm6_chr2L"

    if not dest.exists():
        url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/{suffix}"
        _download_with_retry(url, genome_path / suffix)
        with tarfile.open(genome_path / suffix) as tar:
            dest.mkdir(parents=True, exist_ok=True)
            tar.extractall(path=dest)
        os.remove(genome_path / suffix)

    return dest

@pytest.fixture(scope="session")
def chromsizes(genome_path: Path) -> Path:
    suffix = "chr21.fa.fai"
    dest = genome_path / suffix

    if not (genome_path / "chr21_rename.fa.fai").exists():
        url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/{suffix}"
        _download_with_retry(url, dest)

        # Append dm6 chromsizes
        with open(dest, "ab") as f:
            url2 = (
                "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes"
            )
            for attempt in range(3):
                try:
                    r2 = requests.get(url2, stream=True, timeout=30)
                    r2.raise_for_status()
                    for line in r2.iter_lines():
                        f.write(b"dm6_" + line + b"\n")
                    break
                except requests.RequestException:
                    if attempt < 2:
                        time.sleep(2**attempt)
                        continue
                    raise

    return dest


@pytest.fixture(scope="function")
def gtf(genome_path: Path, assay: str, index: Path) -> Path:
    if "rna" in assay:
        gtf_path = genome_path / "chr21_rna_spikein.gtf"
    else:
        gtf_path = genome_path / "chr21.gtf"

    if not gtf_path.exists():
        url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/{gtf_path.name}"
        _download_with_retry(url, gtf_path)
    return gtf_path


@pytest.fixture(scope="function")
def blacklist(genome_path: Path) -> Path:
    blacklist_path = genome_path / "hg38-blacklist.v2.bed.gz"
    if not blacklist_path.exists():
        url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
        _download_with_retry(url, blacklist_path)
    return blacklist_path


@pytest.fixture(scope="function")
def genome_files(genome_path: Path, assay: str) -> tuple[Path, Path]:
    """Only download genome files if needed for meth/snp assays."""
    fasta = genome_path / "chr21_meth.fa"
    fasta_fai = genome_path / "chr21_meth.fa.fai"

    # Only download if this assay actually needs it
    if assay in ["meth", "snp"]:
        if not fasta.exists():
            url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/chr21_meth.fa"
            _download_with_retry(url, fasta)
        if not fasta_fai.exists():
            url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/chr21_meth.fa.fai"
            _download_with_retry(url, fasta_fai)

    return fasta, fasta_fai


@pytest.fixture(scope="session")
def mcc_files(genome_path: Path) -> dict[str, Path]:
    files: dict[str, Path] = {}
    mcc_viewpoints = genome_path / "mcc_viewpoints.bed"
    if not mcc_viewpoints.exists():
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/seqnado_data/test_viewpoints.bed"
        _download_with_retry(url, mcc_viewpoints)
    files["viewpoints"] = mcc_viewpoints
    return files


# ------------------
# Assay helpers
# ------------------


@pytest.fixture(scope="session")
def cores(pytestconfig: pytest.Config) -> int:
    """Number of cores to use for pipeline executions."""
    return int(pytestconfig.getoption("--cores"))


@pytest.fixture(scope="session")
def selected_assays(pytestconfig: pytest.Config) -> list[str]:
    """List of assays to parametrize pipeline tests with."""
    return _parse_assays(pytestconfig)


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
        tar_path = target_dir.parent / "fastq.tar.gz"
        _download_with_retry(url, tar_path)
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
    # Use getbasetemp() to get shared base temp directory, then create assay-specific subdirectory
    # This ensures all tests for the same assay use the same directory
    base_temp = tmp_path_factory.getbasetemp()
    run_dir = base_temp / assay
    run_dir.mkdir(exist_ok=True)
    return run_dir


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
def run_init(
    star_index: Path,
    bt2_index: Path,
    chromsizes: Path,
    gtf: Path,
    blacklist: Path,
    run_directory: Path,
    assay: str,
    monkeypatch: pytest.MonkeyPatch,
    genome_files: tuple[Path, Path],
    genome_path: Path,
):
    """Run `seqnado init` and materialize genome_config.json with appropriate indices."""
    monkeypatch.setenv("SEQNADO_CONFIG", str(run_directory))
    monkeypatch.setenv("HOME", str(run_directory))

    process = subprocess.run(
        ["seqnado", "init"],
        input="y\n",
        text=True,
        cwd=run_directory,
        capture_output=True,
    )

    genome_config_file = run_directory / ".config" / "seqnado" / "genome_config.json"

    fasta, _ = genome_files
    genome_config_dict = {
        "hg38": {
            "star_index": str(star_index),
            "bt2_index": str(bt2_index),
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
    dest_fastq_screen_config = (
        run_directory / ".config" / "seqnado" / "fastq_screen.conf"
    )
    dest_fastq_screen_config.parent.mkdir(parents=True, exist_ok=True)

    if source_fastq_screen_config.exists():
        # Read template and replace index path using regex for robustness
        with open(source_fastq_screen_config, "r") as f:
            config_content = f.read()
        # Replace DATABASE lines with actual index path
        # Pattern: DATABASE\tName\t/any/path -> DATABASE\tName\t{actual_index}
        config_content = re.sub(
            r"(DATABASE\s+\S+\s+).*", rf"\1{str(bt2_index)}", config_content
        )
        with open(dest_fastq_screen_config, "w") as f:
            f.write(config_content)
    else:
        # Fallback: create minimal config
        with open(dest_fastq_screen_config, "w") as f:
            f.write(f"""# Test fastq_screen.conf file

DATABASE\tTest\t{str(bt2_index)}
""")

    assert process.returncode == 0, f"seqnado init failed: {process.stderr}"
    assert genome_config_file.exists(), "genome_config.json missing after init"


@pytest.fixture(scope="function")
def config_yaml(
    run_directory: Path, assay_type: str, monkeypatch: pytest.MonkeyPatch
) -> Path:
    monkeypatch.setenv("SEQNADO_CONFIG", str(run_directory))
    monkeypatch.setenv("HOME", str(run_directory))

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

    assert result.returncode == 0, (
        f"seqnado config failed: {result.stderr}\n{result.stdout}"
    )
    assert config_file_path.exists(), f"{assay_type} config file not created."
    return config_file_path


@pytest.fixture(scope="function")
def config_yaml_for_testing(config_yaml: Path, assay: str) -> Path:
    import yaml

    with open(config_yaml, "r") as f:
        config = yaml.safe_load(f)

    # set plot_with_plotnado to False
    config["assay_config"]["plot_with_plotnado"] = False

    match assay:
        case "atac":
            config["pileup_method"] = ["deeptools", "homer"]
            config["call_peaks"] = True
            config["peak_calling_method"] = ["lanceotron", "macs", "homer"]
        case "cat":
            config["pileup_method"] = ["deeptools", "bamnado"]
        case "chip":
            config["scale"] = "yes"
            config["library_complexity"] = False
        case "chip-rx":
            config["peak_calling_method"] = ["seacr"]
        case "mcc":
            config["call_peaks"] = False
            config["peak_calling_method"] = ["lanceotron-mcc"]
        case "rna":
            config["pileup_method"] = ["deeptools", "bamnado"]
        case "rna-rx":
            config["pileup_method"] = ["deeptools", "bamnado"]
        case _:
            pass

    with open(config_yaml, "w") as f:
        yaml.dump(config, f)
    return Path(config_yaml)


@pytest.fixture(scope="function")
def seqnado_run_dir(config_yaml_for_testing: Path) -> Path:
    return Path(config_yaml_for_testing).parent


@pytest.fixture(scope="function")
def design(
    seqnado_run_dir: Path, assay_type: str, assay: str, run_directory: Path
) -> Path:
    import pandas as pd

    # Copy FASTQs from run_directory to seqnado_run_dir
    for fq in run_directory.glob("*.fastq.gz"):
        shutil.copy(fq, seqnado_run_dir)

    metadata_file = f"metadata_{assay_type}.csv"
    cmd = [
        "seqnado",
        "design",
        assay_type,
        "-o",
        metadata_file,
        "--no-interactive",
        "--accept-all-defaults",
    ]
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
    # Use environment variable if set, otherwise fall back to home directory
    apptainer_cache_dir = Path(
        os.environ.get("APPTAINER_CACHEDIR")
        or os.environ.get("SINGULARITY_CACHEDIR")
        or (Path.home() / ".apptainer" / "cache")
    )
    multiqc_config = (
        Path(importlib.resources.files(seqnado.data) / "multiqc_config.yaml")
        .absolute()
        .resolve()
    )
    multiqc_config_parent = multiqc_config.parent

    monkeypatch.setenv(
        "APPTAINER_BINDPATH",
        f"{wd}:{wd},{test_data_path}:{test_data_path},{indicies_mount}:{indicies_mount},{tmpdir}:{tmpdir},{multiqc_config_parent}:{multiqc_config_parent}",
    )
    # Set both for compatibility (Apptainer supports both variable names)
    monkeypatch.setenv("APPTAINER_CACHEDIR", str(apptainer_cache_dir), prepend=False)
    monkeypatch.setenv("SINGULARITY_CACHEDIR", str(apptainer_cache_dir), prepend=False)
    monkeypatch.setenv("APPTAINER_TMPDIR", str(tmpdir), prepend=False)


# -------------------------
# Multi-assay test fixtures
# -------------------------


@pytest.fixture(scope="function")
def multi_assay_run_directory(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Create a shared directory for multi-assay tests."""
    return tmp_path_factory.mktemp("multi_assay")


@pytest.fixture(scope="function")
def multi_assay_configs(
    multi_assays: list[str],
    multi_assay_run_directory: Path,
    fastqs: dict[str, list[Path]],
    index: Path,
    chromsizes: Path,
    gtf: Path,
    blacklist: Path,
    genome_path: Path,
    genome_files: tuple[Path, Path],
    monkeypatch: pytest.MonkeyPatch,
):
    """Set up configs and metadata for multiple assays in one directory.

    Returns:
        Dict mapping assay -> {"config": Path, "metadata": Path}
    """
    import yaml

    result = {}
    cwd = Path.cwd()

    try:
        os.chdir(multi_assay_run_directory)

        # Run seqnado init once for the shared directory
        monkeypatch.setenv("SEQNADO_CONFIG", str(multi_assay_run_directory))
        monkeypatch.setenv("HOME", str(multi_assay_run_directory))

        init_process = subprocess.run(
            ["seqnado", "init"],
            input="y\n",
            text=True,
            cwd=multi_assay_run_directory,
            capture_output=True,
        )

        if init_process.returncode != 0:
            print(f"seqnado init STDOUT:\n{init_process.stdout}")
            print(f"seqnado init STDERR:\n{init_process.stderr}")

        assert init_process.returncode == 0, (
            f"seqnado init failed: {init_process.stderr}"
        )

        # Check if .config/seqnado directory was created
        seqnado_dir = multi_assay_run_directory / ".config" / "seqnado"
        if not seqnado_dir.exists():
            print(f"ERROR: .config/seqnado directory not created at {seqnado_dir}")
            print(f"Directory contents: {list(multi_assay_run_directory.iterdir())}")
            raise FileNotFoundError(
                f".config/seqnado directory not created at {seqnado_dir}"
            )

        # Update genome config with test paths - OVERWRITE to exclude dm6 placeholders
        genome_config_path = seqnado_dir / "genome_config.json"

        # Create genome config with ONLY hg38 (overwrite any dm6 from seqnado init)
        genome_config = {
            "hg38": {
                "star_index": str(genome_path / "star_index"),
                "bt2_index": str(index),
                "chromosome_sizes": str(chromsizes),
                "gtf": str(gtf),
                "blacklist": str(blacklist),
                "genes": "",
                "fasta": str(genome_files[0]),
            }
        }

        with open(genome_config_path, "w") as f:
            json.dump(genome_config, f, indent=2)

        # Create fastq_screen.conf file
        fastq_screen_config_path = seqnado_dir / "fastq_screen.conf"
        with open(fastq_screen_config_path, "w") as f:
            f.write(f"""# Test fastq_screen.conf file

DATABASE\tTest\t{index}
""")

        # Copy ALL fastqs for all assays to the shared directory
        all_fastqs = set()
        for assay in multi_assays:
            all_fastqs.update(fastqs.get(assay, []))

        for fq in all_fastqs:
            shutil.copy(fq, multi_assay_run_directory)

        # Create config and metadata for each assay
        for assay in multi_assays:
            # Create config
            config_result = subprocess.run(
                [
                    "seqnado",
                    "config",
                    assay,
                    "--no-interactive",
                    "--no-make-dirs",
                    "--render-options",
                ],
                text=True,
                cwd=multi_assay_run_directory,
                capture_output=True,
            )

            assert config_result.returncode == 0, (
                f"seqnado config {assay} failed: {config_result.stderr}"
            )

            config_path = multi_assay_run_directory / f"config_{assay}.yaml"
            assert config_path.exists(), f"config_{assay}.yaml not created"

            # Modify config for testing
            with open(config_path, "r") as f:
                config = yaml.safe_load(f)

            # Apply test-specific modifications
            config["assay_config"]["create_ucsc_hub"] = False
            config["call_peaks"] = (
                False  # Disable peak calling for tests (test data too small for LanceOtron)
            )

            if assay == "chip":
                config["scale"] = "yes"
                config["library_complexity"] = False
            elif assay == "atac":
                config["pileup_method"] = ["deeptools"]

            with open(config_path, "w") as f:
                yaml.dump(config, f)

            # Create metadata/design
            metadata_file = f"metadata_{assay}.csv"

            # Get list of FASTQs for this assay that were just copied
            assay_fastqs = list(multi_assay_run_directory.glob("*.fastq.gz"))

            design_cmd = [
                "seqnado",
                "design",
                assay,
                "-o",
                metadata_file,
                "--no-interactive",
                "--accept-all-defaults",
                f"fastq/{assay}*.fastq.gz",
            ]
            # Add FASTQ files as arguments
            design_cmd.extend([str(fq) for fq in assay_fastqs])

            design_result = subprocess.run(
                design_cmd,
                cwd=multi_assay_run_directory,
                capture_output=True,
                text=True,
            )
            assert design_result.returncode == 0, (
                f"seqnado design {assay} failed: {design_result.stderr}"
            )

            metadata_path = multi_assay_run_directory / metadata_file
            assert metadata_path.exists(), f"metadata_{assay}.csv not created"

            result[assay] = {
                "config": config_path,
                "metadata": metadata_path,
            }

        yield result

    finally:
        os.chdir(cwd)
