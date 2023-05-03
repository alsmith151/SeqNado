import glob
import os
import shutil
import subprocess
from datetime import datetime
from cookiecutter.main import cookiecutter
import pytest


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
    return glob.glob(os.path.join(path, "rna*.fastq.gz"))


@pytest.fixture(scope="module")
def config_path(data_path):
    return os.path.join(data_path, "config")


@pytest.fixture(scope="module")
def genome_indicies(genome_path):

    indicies = os.path.join(genome_path, "GenomeDir")
    gtf = os.path.join(genome_path, "chr21.gtf")
    fasta = os.path.join(genome_path, "chr21_rename.fa")

    if not os.path.exists(indicies):
        try:
            import requests
            import tarfile

            url = "https://userweb.molbiol.ox.ac.uk/public/asmith/ngs_pipeline/star.tar.gz"
            output = os.path.join(genome_path, "star.tar.gz")
            r = requests.get(url, stream=True)
            with open(output, "wb") as f:
                f.write(r.content)

            tar = tarfile.open(output)
            tar.extractall(path=genome_path)
            tar.close()
            os.remove(output)
            os.rename(os.path.join(genome_path, "GenomeDir"), indicies)

        except Exception as e:
            print(e)
            print("Could not download indicies so generating them")
            os.mkdir(indicies)
            cmd = f"""STAR
                  --runMode genomeGenerate
                  --runThreadN 4
                  --genomeDir {indicies}
                  --genomeFastaFiles {fasta}
                  --sjdbGTFfile {gtf}
                  --sjdbOverhang 100
                  --genomeSAindexNbases 11
                  """
            subprocess.run(cmd.split())

    return indicies


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
    genome_path,
    genome_indicies,
    chromsizes,
    fastqs,
    config_path,
):

    cwd = os.getcwd()
    os.chdir(run_directory)

    cookiecutter(
    f"{package_path}/workflow/config/cookiecutter_config/config_rna/",
    extra_context={
        "genome": "hg19",
        "date": "{% now 'utc', '%Y-%m-%d' %}",
        "project_name": "test",
        "chromosome_sizes": chromsizes,
        "indicies": genome_indicies,
        "remove_blacklist": "yes",
        "blacklist": f"{data_path}/genome/hg19_blacklist.bed",
        "make_ucsc_hub": "yes",
        "UCSC_hub_directory": "test_hub",
        "email": "test",
        "color_by": "samplename",
        "gtf": f"{data_path}/genome/chr21.gtf",
    },
    no_input=True,
    )

    # Move config files and fastq files
    current_date = datetime.now().strftime("%Y-%m-%d")
    os.chdir(f"{current_date}_test")
    for fq in fastqs:
        shutil.copy(fq, ".")


    yield

    os.chdir(cwd)


def test_pipeline_singularity(genome_path, genome_indicies, chromsizes):

    cmd = [
        "seqnado",
        "rna",
        "--cores",
        "4",
        "--configfile",
        "config_rna.yml",
        "--use-singularity",
        "--singularity-args",
        f'" -B {genome_indicies} -B {genome_path} -B {chromsizes} "',
    ]
    completed = subprocess.run(" ".join(cmd), shell=True)
    assert completed.returncode == 0
    assert not os.path.exists("seqnado_error.log")
    assert os.path.exists("seqnado_output/")
