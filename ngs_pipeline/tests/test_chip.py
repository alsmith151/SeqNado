import glob
import os
import shutil
import subprocess

import pytest


@pytest.fixture(scope="module")
def package_path():

    path_file = os.path.abspath(__file__)
    path_test = os.path.dirname(path_file)
    path_package = os.path.dirname(path_test)
    return path_package


@pytest.fixture(scope="module")
def repo_path(package_path):
    return os.path.dirname(package_path)


@pytest.fixture(scope="module")
def pipeline_path(package_path):
    return os.path.join(package_path, "chipseq", "snakefile")


@pytest.fixture(scope="module")
def data_path(package_path):
    return os.path.join(package_path, "data")


@pytest.fixture(scope="module")
def chromsizes(data_path):
    return os.path.join(data_path, "genome", "chr21_rename.fa.fai")


@pytest.fixture(scope="module")
def genome_indicies(data_path):
    indicies = os.path.join(data_path, "genome", "bt2")

    if not os.path.exists(indicies):

        fasta = os.path.join(data_path, "genome", "chr21_rename.fa")
        os.mkdir(indicies)
        cmd = f"bowtie2-build {fasta} {indicies}/chr21 --threads 8"
        subprocess.run(cmd.split())

    return indicies


@pytest.fixture(scope="module")
def run_directory(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data")
    return fn


@pytest.fixture(scope="module", autouse=True)
def set_up(
    run_directory, data_path, package_path, repo_path, genome_indicies, chromsizes
):

    cwd = os.getcwd()
    os.chdir(run_directory)

    # Move config files and fastq files
    for fq in glob.glob(f"{data_path}/test/*CTCF*_*_*.fastq.gz"):
        shutil.copy(fq, ".")

    # Move and replace the config file
    replacements_dict = {
        "GENOME_NAME": "hg19",
        "CHROMOSOME_SIZES_FILE": f"{chromsizes}",
        "INDICES_DIRECTORY": f"{genome_indicies}/chr21",
        "HUB_DIRECTORY_PATH": f"{run_directory}/UCSC_HUB/",
        "MACS_OPTIONS": "-g 46709983",
    }

    with open(f"{repo_path}/config_chip.yml", "r") as r:
        with open("config_chip.yml", "w") as w:
            for line in r:
                for rep_key in replacements_dict:
                    if rep_key in line:
                        line = line.replace(rep_key, replacements_dict[rep_key])

                w.write(line)

    yield

    os.chdir(cwd)


def test_pipeline(pipeline_path):

    cmd = (
        f"ngs-pipeline chip --cores 1 --configfile config_chip.yml"
    )
    completed = subprocess.run(cmd.split())
    assert completed.returncode == 0
