#!/bin/bash

# Set env variables
export APPTAINER_BINDPATH="" # Don't bind any paths as this causes errors

# Build all containers
for container in $(ls *.def);
do
    echo "Building container $container"
    apptainer build --fakeroot $(basename $container .def).sif $container
done