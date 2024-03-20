import glob
import os
import pathlib
import re
import shutil
import subprocess
import tarfile
from datetime import datetime

import pytest
import requests


@pytest.fixture(scope="function", params=["atac", "chip", "chip-rx", "rna", "rna-rx"], autouse=True)
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
    p =  test_data_path / "genome"
    p.mkdir(parents=True, exist_ok=True)
    return p


@pytest.fixture(scope="function")
def genome_indices_path(genome_path, assay) -> pathlib.Path:

    if "rna" not in assay:
        return genome_path / "bt2_chr21_dm6_chr2L"
    else:
        return genome_path / "STAR_chr21_rna_spikein"


@pytest.fixture(scope="function", autouse=True)
def indicies(genome_indices_path, genome_path):

    download_indices = True if not genome_indices_path.exists() else False
    suffix = genome_indices_path.with_suffix(".tar.gz").name
    url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/{suffix}"

    if "bt2" in str(genome_indices_path):
        indicies_path = genome_indices_path / "bt2_chr21_dm6_chr2L"
    else:
        indicies_path = genome_indices_path
    
    if download_indices:

        r = requests.get(url, stream=True)

        tar_index = genome_indices_path.with_suffix(".tar.gz")

        with open(tar_index, "wb") as f:
            f.write(r.content)

        tar = tarfile.open(tar_index)

        if "bt2" in str(genome_indices_path): # These are individual files so need to extract to the indicies folder
            genome_indices_path.mkdir(parents=True, exist_ok=True)
            tar.extractall(path=genome_indices_path, filter="data")
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

    return genome_path / suffix


@pytest.fixture(scope="function")
def gtf(genome_path, assay, indicies):

    if "rna" in assay:
        gtf_path = indicies / "chr21.gtf"
    else:
        gtf_path = genome_path / "chr21.gtf"
    
    if not gtf_path.exists() and "rna" not in assay:
        url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/chr21.gtf"
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
def assay_type(assay):
    return re.sub(r"(.*)\-.*", r"\1", assay)


@pytest.fixture(scope="function")
def snakefile_path(package_path, assay_type):
    return package_path / "workflow" / f"snakefile_{assay_type}"


@pytest.fixture(scope="function")
def fastqs(test_data_path, assay) -> list[pathlib.Path]:
    path = test_data_path / "fastq"

    match assay:
        case "atac":
            files = list(path.glob("atac*.fastq.gz"))
        case "chip":
            files = list(path.glob("CTCF*.fastq.gz"))
            files.append(path / "SINGLE-CTCF_CTCF.fastq.gz")
        case "chip-rx":
            files = list(path.glob("chip-rx*.fastq.gz"))
        case "rna":
            files = list(path.glob("rna_*.fastq.gz"))
        case "rna-rx":
            files = list(path.glob("rna-spikein*.fastq.gz"))

    return files


@pytest.fixture(scope="function")
def run_directory(tmpdir_factory, assay):
    fn = tmpdir_factory.mktemp(assay)
    return fn


@pytest.fixture(scope="function")
def user_inputs(test_data_path, indicies, chromsizes, assay, assay_type, gtf, blacklist):

    defaults = {
        "project_name": "test",
        "genome_name": "hg38",
        "indices": indicies,
        "chromsizes": str(chromsizes),
        "gtf": str(gtf),
        "blacklist": str(blacklist),
        "fastq_screen": "no",
        "remove_blacklist": "yes",
    }

    defaults_atac = {
        "remove_pcr_duplicates": "yes",
        "remove_pcr_duplicates_method": "picard",
        "shift_atac_reads": "yes",
        "make_bigwigs": "yes",
        "pileup_method": "deeptools",
        "make_heatmaps": "yes",
        "scale": "no",
        "call_peaks": "yes",
        "peak_calling_method": "lanceotron",
    }

    defaults_chip = {
        "remove_pcr_duplicates": "yes",
        "remove_pcr_duplicates_method": "picard",
        "spikein": "no",
        "make_bigwigs": "yes",
        "pileup_method": "deeptools",
        "make_heatmaps": "yes",
        "scale": "yes",
        "call_peaks": "yes",
        "peak_calling_method": "lanceotron",
    }

    defaults_chip_rx = {
        "remove_pcr_duplicates": "yes",
        "remove_pcr_duplicates_method": "picard",
        "spikein": "yes",
        "normalisation_method": "orlando",
        "reference_genome": "hg38",
        "spikein_genome": "dm6",
        "make_bigwigs": "yes",
        "pileup_method": "deeptools",
        "scale": "no",
        "make_heatmaps": "yes",
        "call_peaks": "yes",
        "peak_calling_method": "lanceotron",
    }

    defaults_rna = {
        "make_bigwigs": "yes",
        "pileup_method": "deeptools",
        "make_heatmaps": "yes",
        "run_deseq2": "yes",
    }

    defaults_rna_rx = {
        "spikein": "yes",
        "run_deseq2": "yes",
    }

    hub = {
        "make_ucsc_hub": "yes",
        "UCSC_hub_directory": "test_hub",
        "email": "test",
        "color_by": "samplename",
    }

    match assay:
        case "atac":
            return {**defaults, **defaults_atac, **hub}
        case "chip":
            return {**defaults, **defaults_chip, **hub}
        case "chip-rx":
            return {**defaults, **defaults_chip_rx, **hub}
        case "rna":
            return {**defaults, **defaults_rna, **hub}
        case "rna-rx":
            return {**defaults, **defaults_rna_rx, **hub}


@pytest.fixture(scope="function")
def config_yaml(run_directory, user_inputs, assay_type):

    user_inputs = "\n".join([str(v) for v in user_inputs.values()])
    cmd = ["seqnado-config", assay_type]

    # Run the script with subprocess
    process = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=run_directory,
    )

    stdout, stderr = process.communicate(input=user_inputs)

    project_name = "test"
    date = datetime.now().strftime("%Y-%m-%d")
    config_file_path = (
        run_directory / f"{date}_{assay_type}_{project_name}/config_{assay_type}.yml"
    )
    return config_file_path


@pytest.fixture(scope="function")
def seqnado_run_dir(config_yaml):
    return pathlib.Path(config_yaml).parent


@pytest.fixture(scope="function")
def design(seqnado_run_dir, assay_type):
    cmd = ["seqnado-design", assay_type]
    completed = subprocess.run(" ".join(cmd), shell=True, cwd=seqnado_run_dir)
    assert completed.returncode == 0
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


def test_config_generation(config_yaml, assay_type):
    assert os.path.exists(config_yaml), f"{assay_type} config file not created."


def test_pipeline(assay, assay_type, config_yaml, indicies, test_data_path, cores):

    indicies_mount = indicies.parent if not indicies.is_dir() else indicies
    tmpdir = pathlib.Path(os.environ.get("TMPDIR", "/tmp"))

    cmd = [
        "seqnado",
        f"{assay_type}",
        "--cores",
        str(cores),
        "--configfile",
        str(config_yaml),
        "--use-singularity",
        "--singularity-args",
        f'" -B {indicies_mount.resolve()} -B {test_data_path} -B {os.getcwd()} -B {tmpdir}"',
    ]
    completed = subprocess.run(" ".join(cmd), shell=True)
    assert completed.returncode == 0
    assert not os.path.exists("seqnado_error.log")
    assert os.path.exists("seqnado_output/")
