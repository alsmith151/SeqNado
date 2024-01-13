# Installation

The recommended way to install SeqNado is via a combination of conda and pip, this will ensure that all dependencies are installed correctly in a virtual environment. Using the recommended method, the main dependencies of the pipeline are accessed from a singularity container, so it is not necessary to install these manually. If required, conda can be used to install the dependencies locally, but this is not recommended.

## Set-up

### Check conda|mamba is installed

Ideally, ensure that mamba is installed and is working correctly.

```bash
which mamba
```

You should see something like:

```bash
/ceph/project/milne_group/asmith/software/mambaforge/condabin/mamba
```

### Install mamba (if required)

If mamba is not installed, install it using the following command:

!!! Warning
    Ensure that this is installed in your /project directory, not your home directory! These environments can be very large and will fill up your home directory very quickly.

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

## Quick installation

Use the following command to install SeqNado and all dependencies:

```bash
bash <(curl -s https://raw.githubusercontent.com/alsmith151/SeqNado/master/install_seqnado.sh)
```

## Detailed installation for advanced users or troubleshooting

### Create a conda environment

Create a basic conda environment (with pip to install python packages) and activate it.

```bash
mamba create -n seqnado pip
conda activate seqnado
```

### Install the pipeline

#### Install the package from pip (recommended)

```bash
pip install seqnado
```

#### Install from GitHub directly

To install the latest version of the pipeline from GitHub (master branch), use the following command:

```bash
pip install git+https://github.com/alsmith151/SeqNado.git
```

You can also install a specific branch e.g. `develop`:

```bash
pip install git+https://github.com/alsmith151/SeqNado.git@develop
```

#### Install from a local copy of the repository

Clone the repositry and install directly.

```bash
git clone https://github.com/alsmith151/SeqNado.git
cd SeqNado
pip install .
```

### Install the dependencies (if required)

Assuming the conda environment is activated, install the dependencies using the following command:

```bash
wget https://raw.githubusercontent.com/alsmith151/SeqNado/master/environment.yml
mamba env update -f environment.yml
```
