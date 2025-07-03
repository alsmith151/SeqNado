import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tarfile
from datetime import datetime
from unittest.mock import patch

import pytest
import requests
from loguru import logger

from seqnado.config import (
    GenomeConfig,
    WorkflowConfig,
    build_workflow_config,
    create_config,
    load_genome_config,
)

logger.add(sys.stderr, level="DEBUG")


@pytest.fixture(
    scope="function",
    params=["atac", "chip", "chip-rx", "rna", "rna-rx", "snp", "cat", "meth", "mcc"],
    autouse=True,
)
def assay(request):
    return request.param


@pytest.fixture(scope="function")
def repo_path() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[1]


@pytest.fixture(scope="function")
def package_path(repo_path):
    return repo_path / "seqnado"


@pytest.fixture(scope="function")
def test_dir_path(repo_path):
    return repo_path / "tests"


@pytest.fixture(scope="function")
def test_data_path(test_dir_path):
    return test_dir_path / "data"


@pytest.fixture(scope="function")
def workflow_path(package_path):
    return package_path / "workflow"


@pytest.fixture(scope="function")
def config_path(workflow_path):
    return workflow_path / "config"


@pytest.fixture(scope="function")
def genome_path(test_data_path):
    p = test_data_path / "genome"
    p.mkdir(parents=True, exist_ok=True)
    return p


@pytest.fixture(scope="function")
def genome_index_path(genome_path, assay) -> pathlib.Path:
    if "rna" in assay:
        return genome_path / "STAR_chr21_rna_spikein"
    elif "meth" in assay:
        return genome_path / "bt2_chr21_meth"
    else:
        return genome_path / "bt2_chr21_dm6_chr2L"


@pytest.fixture(scope="function")
def index(genome_index_path, genome_path) -> pathlib.Path:
    download_index = True if not genome_index_path.exists() else False
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
        tar_index = genome_index_path.with_suffix(".tar.gz")
        with open(tar_index, "wb") as f:
            f.write(r.content)
        tar = tarfile.open(tar_index)
        if "bt2" in str(genome_index_path):
            genome_index_path.mkdir(parents=True, exist_ok=True)
            tar.extractall(path=genome_index_path, filter="data")
        else:
            tar.extractall(path=genome_path, filter="data")
        tar.close()
        os.remove(tar_index)

    return indicies_path


@pytest.fixture(scope="function")
def chromsizes(genome_path):
    suffix = "chr21.fa.fai"

    if not (genome_path / "chr21_rename.fa.fai").exists():
        url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/{suffix}"
        r = requests.get(url, stream=True)
        with open(genome_path / suffix, "wb") as f:
            f.write(r.content)

            # Add the dm6 chromsizes file to the genome path
            url = (
                "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes"
            )
            r = requests.get(url, stream=True)

            for line in r.iter_lines():
                f.write(b"dm6_" + line + b"\n")

    return genome_path / suffix


@pytest.fixture(scope="function")
def gtf(genome_path, assay, index):
    if "rna" in assay:
        gtf_path = genome_path / "chr21_rna_spikein.gtf"
    else:
        gtf_path = genome_path / "chr21.gtf"

    if not gtf_path.exists():
        url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/{gtf_path.name}"
        r = requests.get(url, stream=True)
        with open(gtf_path, "wb") as f:
            f.write(r.content)

    return gtf_path


@pytest.fixture(scope="function")
def blacklist(genome_path):
    blacklist_path = genome_path / "hg38-blacklist.v2.bed.gz"
    if not blacklist_path.exists():
        url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
        r = requests.get(url, stream=True)
        with open(blacklist_path, "wb") as f:
            f.write(r.content)

    return blacklist_path


@pytest.fixture(scope="function")
def genome_files(genome_path):
    fasta = genome_path / "chr21_meth.fa"
    fasta_fai = genome_path / "chr21_meth.fa.fai"

    if not fasta.exists():
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/chr21_meth.fa"
        r = requests.get(url, stream=True)
        with open(fasta, "wb") as f:
            f.write(r.content)

    if not fasta_fai.exists():
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/chr21_meth.fa.fai"
        r = requests.get(url, stream=True)
        with open(fasta_fai, "wb") as f:
            f.write(r.content)

    return fasta, fasta_fai


@pytest.fixture(scope="function")
def mcc_files(genome_path):
    files = dict()

    url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/seqnado_data/test_viewpoints.bed"
    r = requests.get(url, stream=True)
    mcc_viewpoints = genome_path / "mcc_viewpoints.bed"
    with open(mcc_viewpoints, "wb") as f:
        f.write(r.content)

    files["viewpoints"] = genome_path / "mcc_viewpoints.bed"

    return files


@pytest.fixture(scope="function")
def assay_type(assay):
    return re.sub(r"(.*)\-.*", r"\1", assay)


@pytest.fixture(scope="function")
def snakefile_path(package_path, assay_type):
    return package_path / "workflow" / f"snakefile_{assay_type}"


@pytest.fixture(scope="function", autouse=True)
def fastqs(test_data_path, assay) -> list[pathlib.Path]:
    target_dir = test_data_path / "fastq"

    if not target_dir.exists():
        # Create target directory where fastq files will reside.
        target_dir.mkdir(parents=True, exist_ok=True)

        # Download the tarball.
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/seqnado_data/fastq.tar.gz"
        r = requests.get(url, stream=True)
        tar_path = target_dir.parent / "fastq.tar.gz"
        with open(tar_path, "wb") as f:
            f.write(r.content)

        # Extract the tarball to a temporary directory.
        temp_extract_dir = target_dir.parent / "temp_fastq_extracted"
        temp_extract_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(path=temp_extract_dir)

        # Move all extracted fastq files into the target directory.
        for fastq_file in temp_extract_dir.rglob("*.fastq.gz"):
            shutil.move(str(fastq_file), target_dir / fastq_file.name)

        # Clean up temporary directory and tarball.
        shutil.rmtree(temp_extract_dir)
        tar_path.unlink()

    # Select fastq files based on the assay type.
    match assay:
        case "atac":
            files = list(target_dir.glob("atac*.fastq.gz"))
        case "chip":
            files = list(target_dir.glob("chip-rx_*.fastq.gz"))
        case "chip-rx":
            files = list(target_dir.glob("chip-rx_*.fastq.gz"))
        case "rna":
            files = list(target_dir.glob("rna_*.fastq.gz"))
        case "rna-rx":
            files = list(target_dir.glob("rna-spikein*.fastq.gz"))
        case "snp":
            files = list(target_dir.glob("snp*.fastq.gz"))
        case "cat":
            files = list(target_dir.glob("chip-rx_*.fastq.gz"))
        case "meth":
            files = list(target_dir.glob("meth*.fastq.gz"))
        case "mcc":
            files = list(target_dir.glob("mcc*.fastq.gz"))
        case _:
            raise ValueError(f"Unsupported assay: {assay}")

    return files


@pytest.fixture(scope="function")
def plot_bed(test_data_path):
    return test_data_path / "plotting_coordinates.bed"


@pytest.fixture(scope="function")
def run_directory(tmpdir_factory, assay):
    fn = tmpdir_factory.mktemp(assay)
    return fn


@pytest.fixture(scope="function", autouse=True)
def run_init(
    index, chromsizes, gtf, blacklist, run_directory, assay, monkeypatch, genome_files
):
    """
    Runs seqnado-init before each test inside the test directory.
    Ensures genome_config.json is correctly written.
    """
    monkeypatch.setenv("SEQNADO_CONFIG", str(run_directory))
    monkeypatch.setenv("HOME", str(run_directory))

    # Run seqnado-init
    cmd = ["seqnado-init"]
    process = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=run_directory,
    )
    stdout, stderr = process.communicate(input="y")

    genome_config_file = run_directory / ".config" / "seqnado" / "genome_config.json"
    if "rna" in assay:
        index_dict = {
            "star_index": str(index),
            "bt2_index": "NA",
        }
    else:
        index_dict = {
            "star_index": "NA",
            "bt2_index": str(index),
        }

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
    with open(genome_config_file, "w") as f:
        json.dump(genome_config_dict, f, indent=4)

    assert process.returncode == 0, f"seqnado-init failed with stderr: {stderr}"
    assert genome_config_file.exists(), (
        "genome_config.json is still missing after seqnado-init!"
    )


@pytest.fixture(scope="function")
def user_inputs(test_data_path, assay, assay_type, plot_bed, genome_files, mcc_files):
    fasta, fasta_fai = genome_files
    prompts = {
        "Annotate SNPs?": "no",
        "Bigwig method:": "deeptools",
        "Calculate library complexity?": "yes" if assay in ["cat", "snp"] else "no",
        "Call methylation?": "yes",
        "Call peaks?": "yes",
        "Call SNPs?": "yes" if assay == "snp" else "no",
        "Color by (for UCSC hub):": "samplename",
        "Do you have spikein?": "yes" if "rx" in assay else "no",
        "Duplicates removal method:": "picard",
        "Fastqscreen config path:": "/dummy/fastqscreen.conf",
        "Generate consensus counts from Design merge column? (yes/no)": "yes"
        if assay in ["atac", "chip-rx"]
        else "no",
        "Generate GEO submission files?": "yes"
        if assay in ["chip", "rna", "cat"]
        else "no",
        "Genome?": "hg38",
        "Make Bigwigs?": "yes",
        "Make heatmaps?": "yes" if assay == "atac" else "no",
        "Make UCSC hub?": "yes",
        "Methylation assay:": "taps",
        "Normalisation method:": "orlando",
        "Path to bed file with coordinates for plotting": str(plot_bed)
        if not assay == "snp"
        else "",
        "Path to bed file with genes.": "",
        "Path to reference fasta index:": str(fasta_fai) if assay in ["meth", "snp"] else "dummy_ref.fasta.fai",
        "Path to reference fasta:": "dummy_ref.fasta"
        if assay not in ["meth", "mcc"]
        else str(fasta),
        "Path to SNP database:": "dummy_snp_db",
        "Peak calling method:": "lanceotron",
        "Perform fastqscreen?": "no",
        "Perform plotting?": "yes" if not assay == "snp" else "no",
        "Project name?": "test",
        "Quantification method:": "feature_counts",
        "Reference genome:": "hg38",
        "Remove blacklist regions?": "yes",
        "Remove PCR duplicates?": "no" if assay in ["rna", "rna-rx"] else "yes",
        "Run DESeq2?": "no",
        "Salmon index path:": "dummy_salmon_index",
        "Shift ATAC reads?": "yes",
        "Spikein genome:": "dm6",
        "UCSC hub directory:": "dummy_hub_dir",
        "What is your email address?": "test@example.com",
        "Path to viewpoints file: (default: path/to/viewpoints.bed):": str(
            mcc_files["viewpoints"]
        )
        if assay == "mcc"
        else "",
        "Resolution for MCC cooler files:": "100",
        "Make dataset for ML?": "yes" if assay == "cat" else "no",
        "Use regions BED file?": "no",
        "Binsize for dataset:": "10000" if assay == "cat" else "",
    }

    return prompts


@pytest.fixture(scope="function")
def config_yaml(run_directory, assay_type, monkeypatch, user_inputs):
    import pexpect

    monkeypatch.setenv("SEQNADO_CONFIG", str(run_directory))
    monkeypatch.setenv("HOME", str(run_directory))

    child = pexpect.spawn(
        "seqnado-config", args=[assay_type], encoding="utf-8", cwd=run_directory
    )
    child.logfile = sys.stdout

    input_keys, input_values = list(user_inputs.keys()), list(user_inputs.values())

    while True:
        try:
            index = child.expect_exact(input_keys)
            value = input_values[index]
            child.sendline(value)
        except pexpect.EOF:
            break

    date = datetime.now().strftime("%Y-%m-%d")
    config_file_path = (
        run_directory / f"{date}_{assay_type}_test/config_{assay_type}.yml"
    )

    assert config_file_path.exists(), f"{assay_type} config file not created."
    return config_file_path


@pytest.fixture(scope="function")
def config_yaml_for_testing(config_yaml, assay):
    import yaml

    with open(config_yaml, "r") as f:
        config = yaml.safe_load(f)

    if assay == "chip":
        config["scale"] = "yes"
        config["library_complexity"] = False
    elif assay == "chip-rx":
        config["peak_calling_method"] = "seacr"
    elif assay == "atac":
        config["pileup_method"] = ["deeptools", "homer"]
        config["call_peaks"] = True
        config["peak_calling_method"] = ["lanceotron", "macs", "homer"]
        config["peak_calling_method"] = ["lanceotron", "macs", "homer"]

    with open(config_yaml, "w") as f:
        yaml.dump(config, f)

    return pathlib.Path(config_yaml)


@pytest.fixture(scope="function")
def seqnado_run_dir(config_yaml_for_testing):
    return pathlib.Path(config_yaml_for_testing).parent


@pytest.fixture(scope="function")
def design(seqnado_run_dir, assay_type, assay):
    import pandas as pd

    cmd = ["seqnado-design", assay_type, "--merge"]
    completed = subprocess.run(" ".join(cmd), shell=True, cwd=seqnado_run_dir)
    assert completed.returncode == 0

    if assay == "chip":
        # Add merge column to design file
        df = pd.read_csv(seqnado_run_dir / "design.csv", index_col=0)
        df["merge"] = "MLL-MERGED-TOGETHER"
        df.to_csv(seqnado_run_dir / "design.csv")

    elif assay == "rna-rx":
        # Add deseq2 column to design file

        df = pd.read_csv(seqnado_run_dir / "design.csv", index_col=0)
        df["deseq2"] = df.index.str.split("-").str[-2]
        df.to_csv(seqnado_run_dir / "design.csv")

    return seqnado_run_dir / "design.csv"


@pytest.fixture(scope="function", autouse=True)
def set_up(seqnado_run_dir, fastqs):
    cwd = pathlib.Path(os.getcwd())
    os.chdir(seqnado_run_dir)

    # Move fastqs to run directory
    for fq in fastqs:
        shutil.copy(fq, seqnado_run_dir)

    yield

    os.chdir(cwd)


def test_config(config_yaml, assay_type):
    assert os.path.exists(config_yaml), f"{assay_type} config file not created."


def test_design(design, assay_type):
    assert os.path.exists(design), f"{assay_type} design file not created."


@pytest.fixture(scope="function", autouse=True)
def apptainer_args(index, test_data_path):
    import importlib.resources
    import seqnado.data
    import pathlib
    import os

    indicies_mount = index.parent if not index.is_dir() else index
    tmpdir = pathlib.Path(os.environ.get("TMPDIR", "/tmp") or "/tmp")
    wd = pathlib.Path(os.getcwd()).resolve()
    apptainer_cache_dir = pathlib.Path.home() / ".apptainer"
    multiqc_config = (
        pathlib.Path(importlib.resources.files(seqnado.data) / "multiqc_config.yaml")
        .absolute()
        .resolve()
    )
    multiqc_config_parent = multiqc_config.parent

    os.environ["APPTAINER_BINDPATH"] = (
        f"{wd}:{wd},"
        f"{test_data_path}:{test_data_path},"
        f"{indicies_mount}:{indicies_mount},"
        f"{tmpdir}:{tmpdir},"
        f"{multiqc_config_parent}:{multiqc_config_parent}"
    )

    if not os.environ.get("APPTAINER_CACHEDIR"):
        os.environ["APPTAINER_CACHEDIR"] = str(apptainer_cache_dir)


@pytest.fixture(scope="function")
def test_profile_path(workflow_path):
    return workflow_path / "envs" / "profiles" / "profile_test"


def test_pipeline(
    assay_type,
    config_yaml_for_testing,
    apptainer_args,
    cores,
    design,
    index,
    test_data_path,
    package_path,
    test_profile_path,
):
    subprocess.run(
        [
            "seqnado",
            assay_type,
            "-c",
            str(cores),
            "--configfile",
            str(config_yaml_for_testing),
            "--workflow-profile",
            test_profile_path,
            "--apptainer-prefix",
            "/tmp",
        ],
    )

    assert not os.path.exists("seqnado_error.log")
    assert os.path.exists("seqnado_output/")
    assert os.path.exists("seqnado_output/seqnado_report.html")
