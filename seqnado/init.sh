#!/bin/bash

# This script is used to initialize the environment for the SeqNado project.

if apptainer remote list | grep -q "SylabsCloud"; then
    echo "SylabsCloud is already a configured remote."
else
    echo "Adding SylabsCloud remote to apptainer - This will allow you to download and run the pre-made containers"
    apptainer remote add --no-login SylabsCloud cloud.sylabs.io
    apptainer remote use SylabsCloud
fi

IS_CCB=$(hostname | grep -c "imm-")

if [ -z "$APPTAINER_BINDPATH" ]; then
    if [ $IS_CCB -eq 1 ]; then
        export APPTAINER_BINDPATH="/ceph:/ceph, /project:/project, /databank:/databank"
        echo 'export APPTAINER_BINDPATH="/ceph:/ceph, /project:/project, /databank:/databank"' >> ~/.bashrc
    else
        echo "Please set the APPTAINER_BINDPATH environment variable to bind the necessary directories."
    fi
fi
