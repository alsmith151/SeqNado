# import glob
# import os
# import shutil
# import subprocess
# import pytest
# from datetime import datetime
# from cookiecutter.main import cookiecutter


# @pytest.fixture(scope="module")
# def repo_path():
#     path_file = os.path.abspath(__file__)
#     path_test = os.path.dirname(path_file)
#     path_repo = os.path.dirname(path_test)
#     return path_repo


# @pytest.fixture(scope="module")
# def package_path(repo_path):
#     return os.path.join(repo_path, "seqnado")


# @pytest.fixture(scope="module")
# def pipeline_path(package_path):
#     return os.path.join(package_path, "workflow", "snakefile_snp")


# @pytest.fixture(scope="module")
# def test_dir_path(repo_path):
#     return os.path.join(repo_path, "tests")


# @pytest.fixture(scope="module")
# def data_path(test_dir_path):
#     return os.path.join(test_dir_path, "data")


# @pytest.fixture(scope="module")
# def genome_path(data_path):
#     return os.path.join(data_path, "genome")


# @pytest.fixture(scope="module")
# def chromsizes(genome_path):
#     return os.path.join(genome_path, "chr21_rename.fa.fai")


# @pytest.fixture(scope="module")
# def fastqs(data_path):
#     path = os.path.join(data_path, "fastq")
#     return glob.glob(os.path.join(path, "snp*.fastq.gz"))


# @pytest.fixture(scope="module")
# def config_path(data_path):
#     return os.path.join(data_path, "config")


# @pytest.fixture(scope="module")
# def genome_indicies(genome_path):

#     indicies = os.path.join(genome_path, "bt2")

#     if not os.path.exists(indicies):
#         try:
#             import requests
#             import tarfile

#             url = (
#                 "https://userweb.molbiol.ox.ac.uk/public/asmith/ngs_pipeline/bt2.tar.gz"
#             )
#             output = os.path.join(genome_path, "bt2.tar.gz")
#             r = requests.get(url, stream=True)
#             with open(output, "wb") as f:
#                 f.write(r.content)

#             tar = tarfile.open(output)
#             tar.extractall(path=genome_path)
#             tar.close()
#             os.remove(output)
#             os.rename(genome_path + "/bt2", indicies)

#         except Exception as e:
#             print(e)
#             print("Could not download indicies so generating them")
#             os.mkdir(indicies)
#             cmd = f"bowtie2-build {os.path.join(genome_path,'chr21_rename.fa')} {indicies}/bt2 --threads 8"
#             subprocess.run(cmd.split())

#     return os.path.join(indicies, "chr21")


# @pytest.fixture(scope="module")
# def run_directory(tmpdir_factory):
#     fn = tmpdir_factory.mktemp("data")
#     return fn


# @pytest.fixture(scope="module", autouse=True)
# def set_up(
#     run_directory,
#     data_path,
#     package_path,
#     repo_path,
#     genome_indicies,
#     chromsizes,
#     fastqs,
#     config_path,
# ):

#     cwd = os.getcwd()
#     os.chdir(run_directory)

#     cookiecutter(
#     f"{package_path}/workflow/config/cookiecutter_config/config_snp/",
#     extra_context={
#         "genome": "hg19",
#         "date": "{% now 'utc', '%Y-%m-%d' %}",
#         "project_name": "test",
#         "chromosome_sizes": chromsizes,
#         "indicies": genome_indicies,
#         "design": "design.csv",
#         "read_type": "paired",
#         "remove_pcr_duplicates_method": "picard",
#         "remove_blacklist": "yes",
#         "blacklist": f"{data_path}/genome/hg19-blacklist.v2.chr21.bed.gz",
#         "fasta": f"{data_path}/genome/chr21_rename.fa",
#         "fasta_index": f"{data_path}/genome/chr21_rename.fa.fai",
#         "call_snps": "yes",
#         "annotate_snps": "no",
#     },
#     no_input=True,
#     )

#     # Move config files and fastq files
#     current_date = datetime.now().strftime("%Y-%m-%d")
#     os.chdir(f"{current_date}_test")
#     for fq in fastqs:
#         shutil.copy(fq, ".")


#     yield

#     os.chdir(cwd)


# def test_pipeline_singularity(genome_path):
#     indicies_dir = os.path.join(genome_path, "bt2")

#     cmd = [
#         "seqnado",
#         "snp",
#         "--cores",
#         "4",
#         "--configfile",
#         "config_snp.yml",
#         "--use-singularity",
#         "--singularity-args",
#         f'" -B {indicies_dir} -B {genome_path}"',
#     ]
#     completed = subprocess.run(" ".join(cmd), shell=True)
#     assert completed.returncode == 0
#     assert not os.path.exists("seqnado_error.log")
#     assert os.path.exists("seqnado_output/")
