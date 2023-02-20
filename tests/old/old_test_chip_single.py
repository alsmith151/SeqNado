import os
import pytest
import wget
import zipfile
import glob
import shutil
import subprocess

path_file = os.path.abspath(__file__)
dir_tests = os.path.dirname(path_file)
dir_tests_run = os.path.join(dir_tests, "pipeline_chip_single")
dir_package = os.path.dirname(dir_tests)
dir_repo = os.path.dirname(dir_package)
dir_data_test = os.path.join(dir_package, "data", "test")
dir_data_genome = os.path.join(dir_package, "data", "genome")


@pytest.fixture(scope="module", autouse=True)
def set_up():

    path_root = os.getcwd()

    # Make indicies
    indicies = os.path.join(dir_data_genome, "bt2")
    fasta = os.path.join(dir_data_genome, "chr21_rename.fa")
    if not os.path.exists(indicies):
        os.mkdir(indicies)
        cmd = f"bowtie2-build {fasta} {indicies}/chr21 --threads 8"
        subprocess.run(cmd.split())

    # Change wd
    if not os.path.exists(dir_tests_run):
        os.mkdir(dir_tests_run)

    os.chdir(dir_tests_run)

    # Move config files and fastq files
    for fq in glob.glob(f"{dir_data_test}/SINGLE*.fastq.gz"):
        shutil.copy(fq, ".")

    # Move and replace the config file
    replacements_dict = {
        "GENOME_NAME": "hg19",
        "CHROM_SIZES": f"{fasta}.fai",
        "PATH_TO_ALIGNER_INDICIES": f"{indicies}/chr21",
        "HUB_DIR": dir_tests_run,
        "BAM_COVERAGE_OPTIONS": "",
        "MACS_OPTIONS": "-g 46709983",
    }

    with open(f"{dir_repo}/config_chip.yml", "r") as r:
        with open("config_chip.yml", "w") as w:
            for line in r:
                for rep_key in replacements_dict:
                    if rep_key in line:
                        line = line.replace(rep_key, replacements_dict[rep_key])

                w.write(line)

    yield

    os.chdir(path_root)
    shutil.rmtree(dir_tests_run)


def test_pipeline():

    cmd = "seqnado chip make --local -p 4"
    completed = subprocess.run(cmd.split())
    assert completed.returncode == 0
