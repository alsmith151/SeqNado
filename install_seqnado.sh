#!/bin/bash

if ! command -v mamba &> /dev/null; then
  echo "mamba not found. Installing mambaforge..."
  # Install mambaforge using the official installer
  wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -O mambaforge.sh
  bash mambaforge.sh -b -p $HOME/mambaforge
  rm mambaforge.sh
  # Add mambaforge to the PATH
  export PATH="$HOME/mambaforge/bin:$PATH"
fi

# Create a new conda environment
echo "Creating a new conda environment..."
mamba create -n seqnado "python>3.10" pip -y

# Activate the conda environment
echo "Activating the conda environment..."
conda activate seqnado

# Install the seqnado package using pip
echo "Installing seqnado..."
pip install seqnado

# Deactivate the conda environment
echo "Deactivating the conda environment..."
conda deactivate

echo "Done!"
