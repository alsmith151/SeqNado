import glob
import os
import pathlib
import pytest
import shutil
import subprocess
from datetime import datetime


@pytest.fixture(scope="module")
def repo_path():
    path_file = os.path.abspath(__file__)
    path_test = os.path.dirname(path_file)
    path_repo = os.path.dirname(path_test)
    return path_repo


@pytest.fixture(scope="module")
def package_path(repo_path):
    return os.path.join(repo_path, "seqnado")


@pytest.fixture(scope="module")
def pipeline_path(package_path):
    return os.path.join(package_path, "workflow", "snakefile_atac")


@pytest.fixture(scope="module")
def test_dir_path(repo_path):
    return os.path.join(repo_path, "tests")


@pytest.fixture(scope="module")
def data_path(test_dir_path):
    return os.path.join(test_dir_path, "data")


@pytest.fixture(scope="module")
def genome_path(data_path):
    return os.path.join(data_path, "genome")


@pytest.fixture(scope="module")
def chromsizes(genome_path):
    return os.path.join(genome_path, "chr21_rename.fa.fai")


@pytest.fixture(scope="module")
def fastqs(data_path):
    path = os.path.join(data_path, "fastq")
    return glob.glob(os.path.join(path, "atac*.fastq.gz"))


@pytest.fixture(scope="module")
def config_path(data_path):
    return os.path.join(data_path, "config")


@pytest.fixture(scope="module")
def genome_indicies(genome_path):
    indicies = os.path.join(genome_path, "bt2")

    if not os.path.exists(indicies):
        try:
            import requests
            import tarfile

            url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/bt2.tar.gz"
            output = os.path.join(genome_path, "bt2.tar.gz")
            r = requests.get(url, stream=True)
            with open(output, "wb") as f:
                f.write(r.content)

            tar = tarfile.open(output)
            tar.extractall(path=genome_path)
            tar.close()
            os.remove(output)
            os.rename(genome_path + "/bt2", indicies)

        except Exception as e:
            print(e)
            print("Could not download indicies so generating them")
            os.mkdir(indicies)
            cmd = f"bowtie2-build {os.path.join(genome_path,'chr21_rename.fa')} {indicies}/bt2 --threads 8"
            subprocess.run(cmd.split())

    return os.path.join(indicies, "chr21")


@pytest.fixture(scope="module")
def run_directory(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data")
    return fn


@pytest.fixture(scope="module")
def user_inputs(
    data_path,
    genome_indicies,
    chromsizes,
):
    return {
        "project_name": "test",
        "genome_name": "hg19",
        "index": genome_indicies,
        "chromsizes": chromsizes,
        "gtf": f"{data_path}/genome/chr21.gtf",
        "blacklist": f"{data_path}/genome/hg19-blacklist.v2.chr21.bed.gz",
        "read_type": "paired",
        "remove_blacklist": "yes",
        "remove_pcr_duplicates": "yes",
        "remove_pcr_duplicates_method": "picard",
        "shift_atac_reads": "yes",
        "split_fastq": "no",
        "make_bigwigs": "yes",
        "pileup_method": "deeptools",
        "make_heatmaps": "yes",
        "call_peaks": "yes",
        "peak_calling_method": "lanceotron",
        "make_ucsc_hub": "yes",
        "UCSC_hub_directory": "test_hub",
        "email": "test",
        "color_by": "samplename",
    }

@pytest.fixture(scope="module")
def test_seqnado_config_creation(
    run_directory,
    user_inputs
    ):
    temp_dir = pathlib.Path(run_directory)
    date = datetime.now().strftime("%Y-%m-%d")
    config_file_path = temp_dir / f"{date}_atac_test/config_atac.yml"
    user_inputs = "\n".join(user_inputs.values())

    cmd = [
        "seqnado-config", 
        "atac"
    ]        

    # Run the script with subprocess
    process = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=temp_dir
    )

    stdout, stderr = process.communicate(input=user_inputs)

    # Assert that the config file was created
    assert os.path.exists(config_file_path), "Config file not created."

@pytest.fixture(scope="module", autouse=True)
def set_up(
    run_directory,
    fastqs,
    user_inputs,
    test_seqnado_config_creation,
):
    cwd = os.getcwd()
    os.chdir(run_directory)

    # Move config files and fastq files
    current_date = datetime.now().strftime("%Y-%m-%d")
    os.chdir(f"{current_date}_atac_test")
    for fq in fastqs:
        shutil.copy(fq, ".")

    yield

    os.chdir(cwd)


def test_pipeline_singularity(genome_path, cores):
    indicies_dir = os.path.join(genome_path, "bt2")

    cmd = [
        "seqnado",
        "atac",
        "--cores",
        str(cores),
        "--configfile",
        "config_atac.yml",
        "--use-singularity",
        "--singularity-args",
        f'" -B {indicies_dir} -B {genome_path}"',
    ]
    completed = subprocess.run(" ".join(cmd), shell=True)
    assert completed.returncode == 0
    assert not os.path.exists("seqnado_error.log")
    assert os.path.exists("seqnado_output/")
