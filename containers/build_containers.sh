#!/bin/bash

# Set env variables
export APPTAINER_BINDPATH="" # Don't bind any paths as this causes errors

# Build all containers
# Note: --ignore-fakeroot-command may be needed
for container in $(ls *.def);
do
    echo "Building container $container"
    apptainer build $(basename $container .def).sif $container
done

## To upload:
#apptainer sign IMG.sif
#apptainer push IMG.sif library://asmith151/seqnado/XXXXX:TAG
