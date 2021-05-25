from setuptools import setup, find_packages

setup(
name='ngs_pipeline',
version='0.0.1',
author='asmith',
author_email='alastair.smith@ndcls.ox.ac.uk',
packages=find_packages(),
entry_points={'console_scripts': ['ngs-pipeline = ngs_pipeline:cli']},
include_package_data=True,
url='https://github.com/alsmith151/ngs_pipeline',
license='LICENSE',
description='Chipseq pipeline based on cgat-core (ruffus)',
long_description=open('README.md').read(),
long_description_content_type="text/markdown",
python_requires='>=3.8',
install_requires=['paramiko>=2.7.1',
                  'sqlalchemy>=1.3.18',
                  'cgatcore>=0.6.7',
                  'apsw',
                  'ruffus',
                  'drmaa',
                  'click',
                  'trackhub',
                  'seaborn',
                  'deeptools',
                  'wget',
                  'pybedtools']
)