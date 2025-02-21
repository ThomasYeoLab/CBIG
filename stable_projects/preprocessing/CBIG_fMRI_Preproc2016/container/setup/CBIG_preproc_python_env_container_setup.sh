#!/bin/bash

# This script sets up CBIG's Python environment for containers.
# Compared to the original version ($CBIG_CODE_DIR/setup/python_env_setup/CBIG_python_env_generic_setup.sh):
# 1. This script removes redundant codes for containers (e.g., checking the existence of CBIG_py3).
# 2. This script creates the environment with the CBIG_python3_container_env.yml,
#    which removes packages not used in the preprocessing pipeline such as PyTorch.
#
# Written by Tian Fang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

echo "Executing this script will install Python packages used by CBIG."

# store current directory
WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INSTALLER="${WORKDIR}/Miniconda3-latest-Linux-x86_64.sh"

# Set LC_ALL=C to avoid possible errors when install packages via pip
# Reference: https://stackoverflow.com/questions/36394101/pip-install-locale-error-unsupported-locale-setting
export LC_ALL=C

# default installation directory
INSTALLATION_DIR=/app/miniconda
echo "-- Downloading Environment manager (conda) --"
echo "After that, its path will be added into your .bashrc file."
echo "conda's path: $INSTALLATION_DIR"

# download and install miniconda
if [ ! -f $INSTALLER ]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $INSTALLER
    bash $INSTALLER -b -p $INSTALLATION_DIR
else
    bash $INSTALLER -b -p $INSTALLATION_DIR
fi
echo "Existing ~/.bashrc willl be backed up as ~/.bashrc.python_bak"
rsync -az ~/.bashrc{,.python_bak}
echo "export PATH=$INSTALLATION_DIR/bin:"'$PATH' >>~/.bashrc
source ~/.bashrc

# create new CBIG_py3
echo "-- Creating CBIG_py3 environment --"
conda env create --file $WORKDIR/CBIG_python3_container_env.yml

conda_failed=$?

if [ $conda_failed -eq 0 ]; then
    echo "Success!"
    rm -f $INSTALLER
else
    echo "Failed to install"
    rm -r $INSTALLATION_DIR
    rm ~/.bashrc
    mv ~/.bashrc{.python_bak,}
    echo "Reverted back to original state"

fi
