[Back to Index](index.md)

# Installation

SeqNado can be installed using either `mamba` or `pip`. Follow the steps below to set up the package:

## Prerequisites
Before installing SeqNado, ensure that you have the following tools installed on your system:

- **Mamba/Conda**: Required for managing environments and dependencies. Install Conda from [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/), and then install Mamba:
  ```bash
  conda install -n base -c conda-forge mamba
  ```

- **Pip**: Comes pre-installed with Python. Verify installation:
  ```bash
  pip --version
  ```

If these tools are not installed, follow the links provided to set them up before proceeding with SeqNado installation.

## Install from Bioconda
SeqNado is available on Bioconda. To install:
```bash
mamba install -c bioconda seqnado
```

## Install from PyPI
SeqNado is also available on PyPI. To install:
```bash
pip install seqnado
```

## Next Steps

Once SeqNado is installed continue to: 

[Initialize](initialisation.md)