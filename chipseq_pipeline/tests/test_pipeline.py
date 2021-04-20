import os
import pytest
import wget
import zipfile
import glob
import shutil
import subprocess

path_file = os.path.abspath(__file__)
dir_tests = os.path.dirname(path_file)
dir_tests_run = os.path.join(dir_tests, 'pipeline')
dir_package = os.path.dirname(dir_tests)
dir_repo = os.path.dirname(dir_package)
dir_data = os.path.join(dir_package, 'data', 'test')

bt2_ind_url = 'https://genome-idx.s3.amazonaws.com/bt/grch38_1kgmaj.zip'
bt2_ind_zip = os.path.join(dir_data, 'bt2_indicies.zip')
bt2_indicies = glob.glob(f'{dir_data}/*.bt2')


@pytest.fixture(scope='module', autouse=True)
def set_up():

    path_root = os.getcwd()

    # Get indicies
    if not os.path.exists(bt2_ind_zip) and not len(bt2_indicies) > 0:
        wget.download(bt2_ind_url, out=bt2_ind_zip)
    
    if not len(bt2_indicies) > 0:
        with zipfile.ZipFile(bt2_ind_zip, 'r') as zip_ref:
            zip_ref.extractall(dir_data)
    
    # Change wd
    os.chdir(dir_tests_run)
    
    # Move config files and fastq files
    for fq in glob.glob(f'{dir_data}/*.fastq.gz'):
        shutil.copy(fq, '.')
    

    # Move and replace the config file
    replacements_dict = {'GENOME_NAME': 'hg38',
                         'PATH_TO_ALIGNER_INDICIES': f'{dir_data}/grch38_1kgmaj',
                         'HUB_DIR': dir_tests_run}


    with open(f'{dir_repo}/config.yml', 'r') as r:
        with open('config.yml', 'w') as w:
            for line in r:
                for rep_key in replacements_dict:
                    if rep_key in line:
                        line = line.replace(rep_key, replacements_dict[rep_key])
            
                w.write(line)
    
    
    yield

    os.chdir(path_root)
    shutil.rmtree(dir_tests_run)

    
    




def test_pipeline():
    
    cmd = 'chipseq_pipeline make --local -p 4'
    completed = subprocess.run(cmd, shell=True)
    
    assert completed.returncode == 0
    


