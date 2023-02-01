import glob
import os
import shutil
import subprocess

import pytest


@pytest.fixture(scope="module")
def repo_path():
    path_file = os.path.abspath(__file__)
    path_test = os.path.dirname(path_file)
    path_repo = os.path.dirname(path_test)
    return path_repo


@pytest.fixture(scope="module")
def package_path(repo_path):
    return os.path.join(repo_path, "ngs_pipeline")


@pytest.fixture(scope="module")
def pipeline_path(package_path):
    return os.path.join(package_path, "chipseq", "snakefile")


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
def config_path(repo_path):
    return os.path.join(repo_path, "config")


@pytest.fixture(scope="module")
def genome_indicies(genome_path):

    indicies = os.path.join(genome_path, "bt2")

    if not os.path.exists(indicies):
        try:
            import requests
            import tarfile

            url = (
                "https://userweb.molbiol.ox.ac.uk/public/asmith/ngs_pipeline/bt2.tar.gz"
            )
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


@pytest.fixture(scope="module", autouse=True)
def set_up(
    run_directory,
    data_path,
    package_path,
    repo_path,
    genome_indicies,
    chromsizes,
    fastqs,
    config_path,
):

    cwd = os.getcwd()
    os.chdir(run_directory)

    # Move config files and fastq files
    for fq in fastqs:
        shutil.copy(fq, ".")

    # Move and replace the config file
    replacements_dict = {
        "GENOME_NAME": "hg19",
        "CHROMOSOME_SIZES_FILE": chromsizes,
        "INDICES_DIRECTORY": genome_indicies,
        "HUB_DIRECTORY_PATH": f"{run_directory}/UCSC_HUB/",
        "MACS_OPTIONS": "-g 46709983",
    }

    with open(f"{config_path}/config_atac.yml", "r") as r:
        with open("config_atac.yml", "w") as w:
            for line in r:
                for rep_key in replacements_dict:
                    if rep_key in line:
                        line = line.replace(rep_key, replacements_dict[rep_key])

                w.write(line)

    yield

    os.chdir(cwd)


# def test_pipeline_conda():

#     cmd = f"ngs-pipeline atac --cores 4 --configfile config_atac.yml"
#     completed = subprocess.run(cmd.split())
#     assert completed.returncode == 0


def test_pipeline_singularity(genome_path):
    indicies_dir = os.path.join(genome_path, "bt2")

    cmd = [
        "ngs-pipeline",
        "atac",
        "--cores",
        "4",
        "--configfile",
        "config_atac.yml",
        "--use-singularity",
        "--singularity-args",
        f'" -B {indicies_dir} -B {genome_path}"',
    ]
    completed = subprocess.run(" ".join(cmd), shell=True)
    assert completed.returncode == 0
