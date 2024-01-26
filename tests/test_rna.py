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
    return os.path.join(package_path, "workflow", "snakefile_rna")


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
def genome_indices(genome_path):
    indices = os.path.join(genome_path, "GenomeDir")
    gtf = os.path.join(genome_path, "chr21.gtf")
    fasta = os.path.join(genome_path, "chr21_rename.fa")

    if not os.path.exists(indices):
        try:
            import requests
            import tarfile

            url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/star.tar.gz"
            output = os.path.join(genome_path, "star.tar.gz")
            r = requests.get(url, stream=True)
            with open(output, "wb") as f:
                f.write(r.content)

            tar = tarfile.open(output)
            tar.extractall(path=genome_path)
            tar.close()
            os.remove(output)
            os.rename(os.path.join(genome_path, "GenomeDir"), indices)

        except Exception as e:
            print(e)
            print("Could not download indices so generating them")
            os.mkdir(indices)
            cmd = f"""STAR
                  --runMode genomeGenerate
                  --runThreadN 4
                  --genomeDir {indices}
                  --genomeFastaFiles {fasta}
                  --sjdbGTFfile {gtf}
                  --sjdbOverhang 100
                  --genomeSAindexNbases 11
                  """
            subprocess.run(cmd.split())

    return indices


@pytest.fixture(scope="module")
def run_directory(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data")
    return fn



@pytest.fixture(scope="module")
def user_inputs(
    data_path,
    genome_indices,
    chromsizes,
):
    return {
        "project_name": "test",
        "genome_name": "hg19",
        "indices": genome_indices,
        "chromsizes": chromsizes,
        "gtf": f"{data_path}/genome/chr21.gtf",
        "blacklist": f"{data_path}/genome/hg19-blacklist.v2.chr21.bed.gz",
        "remove_blacklist": "yes",
        "remove_pcr_duplicates": "no",
        "make_bigwigs": "yes",
        "pileup_method": "deeptools",
        "make_heatmaps": "yes",
        "run_DESeq2": "no",
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
    config_file_path = temp_dir / f"{date}_rna_test/config_rna.yml"
    user_inputs = "\n".join(user_inputs.values())

    cmd = [
        "seqnado-config", 
        "rna"
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
    test_seqnado_config_creation
):
    cwd = os.getcwd()
    os.chdir(run_directory)

    # Move config files and fastq files
    current_date = datetime.now().strftime("%Y-%m-%d")
    os.chdir(f"{current_date}_rna_test")
    for fq in fastqs:
        shutil.copy(fq, ".")

    yield

    os.chdir(cwd)


def test_pipeline_singularity(genome_path, genome_indices, chromsizes, cores):
    cmd = [
        "seqnado",
        "rna",
        "--cores",
        str(cores),
        "--configfile",
        "config_rna.yml",
        "--use-singularity",
        "--singularity-args",
        f'" -B {genome_indices} -B {genome_path} -B {chromsizes} "',
    ]
    completed = subprocess.run(" ".join(cmd), shell=True)
    assert completed.returncode == 0
    assert not os.path.exists("seqnado_error.log")
    assert os.path.exists("seqnado_output/")
