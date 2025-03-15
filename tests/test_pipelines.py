import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tarfile
from datetime import datetime

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
    params=["atac", "chip", "chip-rx", "rna", "rna-rx", "snp", "cat", "meth"],
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
    if "rna" or "meth" not in assay:
        return genome_path / "bt2_chr21_dm6_chr2L"
    elif "meth" in assay:
        return genome_path / "hg38_taps_spikein"
    else:
        return genome_path / "STAR_chr21_rna_spikein"


@pytest.fixture(scope="function")
def index(genome_index_path, genome_path) -> pathlib.Path:
    download_index = True if not genome_index_path.exists() else False
    suffix = genome_index_path.with_suffix(".tar.gz").name
    url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/{suffix}"

    if "bt2" in str(genome_index_path):
        indicies_path = genome_index_path / "bt2_chr21_dm6_chr2L"
    elif "taps" in str(genome_index_path):
        indicies_path = genome_index_path / "bt2_chr21_meth"
    else:
        indicies_path = genome_index_path
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
def meth_genome(genome_path):
    meth_genome_path = genome_path / "chr21_meth"
    meth_genome_path.mkdir(parents=True, exist_ok=True)

    fasta_file = meth_genome_path / "chr21_meth.fa"
    fai_file = meth_genome_path / "chr21_meth.fa.fai"

    if not fasta_file.exists():
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/chr21_meth.fa"
        r = requests.get(url, stream=True)
        if r.status_code == 200:
            with open(fasta_file, "wb") as f:
                f.write(r.content)
        else:
            pytest.fail(f"Failed to download {url}: {r.status_code}")

    if not fai_file.exists():
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/chr21_meth.fa.fai"
        r = requests.get(url, stream=True)
        if r.status_code == 200:
            with open(fai_file, "wb") as f:
                f.write(r.content)
        else:
            pytest.fail(f"Failed to download {url}: {r.status_code}")

    return meth_genome_path


@pytest.fixture(scope="function")
def assay_type(assay):
    return re.sub(r"(.*)\-.*", r"\1", assay)


@pytest.fixture(scope="function")
def snakefile_path(package_path, assay_type):
    return package_path / "workflow" / f"snakefile_{assay_type}"


@pytest.fixture(scope="function")
def fastqs(test_data_path, assay) -> list[pathlib.Path]:
    path = test_data_path / "fastq"

    if not path.exists():
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/fastq.tar.gz"
        r = requests.get(url, stream=True)

        tar_path = path.with_suffix(".tar.gz")

        with open(tar_path, "wb") as f:
            f.write(r.content)

        with tarfile.open(tar_path) as tar:
            tar.extractall(path=path.parent, filter="data")

    match assay:
        case "atac":
            files = list(path.glob("atac*.fastq.gz"))
        case "chip":
            files = list(path.glob("chip-rx_*.fastq.gz"))
        case "chip-rx":
            files = list(path.glob("chip-rx_*.fastq.gz"))
        case "rna":
            files = list(path.glob("rna_*.fastq.gz"))
        case "rna-rx":
            files = list(path.glob("rna-spikein*.fastq.gz"))
        case "snp":
            files = list(path.glob("snp*.fastq.gz"))
        case "cat":
            files = list(path.glob("chip-rx_*.fastq.gz"))
        case "meth":
            files = list(path.glob("meth*.fastq.gz"))

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
    index, chromsizes, gtf, blacklist, run_directory, assay, meth_genome, monkeypatch
):
    monkeypatch.setenv("SEQNADO_CONFIG", str(run_directory))
    monkeypatch.setenv("HOME", str(run_directory))

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

    index_dict = {
        "star_index": "NA",
        "bt2_index": str(meth_genome) if "meth" in assay else str(index),
    }

    genome_config_dict = {
        "hg38": {
            "star_index": index_dict["star_index"],
            "bt2_index": index_dict["bt2_index"],
            "chromosome_sizes": str(chromsizes),
            "gtf": str(gtf),
            "blacklist": str(blacklist),
        }
    }

    with open(genome_config_file, "w", encoding="utf-8") as f:
        json.dump(genome_config_dict, f, indent=4)

    assert process.returncode == 0, f"seqnado-init failed with stderr: {stderr}"
    assert genome_config_file.exists(), "genome_config.json missing after seqnado-init!"


@pytest.fixture(scope="function")
def user_inputs(test_data_path, assay, assay_type, plot_bed):
    if assay == "meth":
        meth_fasta = str(pathlib.Path(test_data_path / "genome" / "chr21_meth.fa"))
        meth_fai = str(pathlib.Path(test_data_path / "genome" / "chr21_meth.fa.fai"))

    prompts = {
        "Bigwig method:": "deeptools",
        "Calculate library complexity?": "yes" if assay == "atac" else "no",
        "Call methylation?": "yes",
        "Call peaks?": "yes",
        "Call SNPs?": "no",
        "Color by (for UCSC hub):": "samplename",
        "Do you have spikein? (yes/no)": "yes" if "rx" in assay else "no",
        "Duplicates removal method:": "picard",
        "Fastqscreen config path:": "/dummy/fastqscreen.conf",
        "Generate consensus counts from Design merge column? (yes/no)": "yes"
        if assay == ["chip-rx", "atac"]
        else "no",
        "Generate GEO submission files?": "yes" if assay in ["chip", "rna"] else "no",
        "Genome?": "hg38",
        "Make Bigwigs?": "yes",
        "Make heatmaps?": "yes" if assay == "atac" else "no",
        "Make UCSC hub?": "yes",
        "Methylation assay": "taps",
        "Normalisation method:": "orlando",
        "Path to bed file with coordinates for plotting": str(plot_bed)
        if not assay == "snp"
        else "",
        "Path to bed file with genes.": "",
        "Path to reference fasta:": "dummy_ref.fasta"
        if not assay == "meth"
        else meth_fasta,
        "Path to reference fasta index:": "dummy_ref.fasta.fai"
        if not assay == "meth"
        else meth_fai,
        "Path to SNP database:": "dummy_snp_db",
        "Peak calling method:": "lanceotron",
        "Perform fastqscreen?": "no",
        "Perform plotting?": "yes" if not assay == "snp" else "no",
        "Project name?": "test",
        "Quantification method:": "feature_counts",  # default RNA response
        "Reference genome:": "hg38",
        "Remove blacklist regions?": "yes",
        "Remove PCR duplicates?": "yes",
        "Run DESeq2?": "no",
        "Salmon index path:": "dummy_salmon_index",
        "Shift ATAC reads?": "yes",
        "SNP caller:": "bcftools",
        "Spikein genome:": "dm6",
        "UCSC hub directory:": "dummy_hub_dir",
        "What is your email address?": "test@example.com",
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
    indicies_mount = index.parent if not index.is_dir() else index
    tmpdir = pathlib.Path(os.environ.get("TMPDIR", "/tmp") or "/tmp")
    wd = pathlib.Path(os.getcwd()).resolve()
    apptainer_cache_dir = pathlib.Path.home() / ".apptainer"
    os.environ["APPTAINER_BINDPATH"] = (
        f"{wd}:{wd},"
        f"{test_data_path}:{test_data_path},"
        f"{indicies_mount}:{indicies_mount},"
        f"{tmpdir}:{tmpdir}"
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
    assert os.path.exists("seqnado_output/qc/full_qc_report.html")
