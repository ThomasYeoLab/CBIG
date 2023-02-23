#! /bin/bash

# Written by Wei Hou Tan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

echo "Executing this script will install Python packages used by CBIG."
read -p "Are you sure? (y/n) " answer

if echo "$answer" | grep -iq "^y"; then
    # store current directory
    WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    INSTALLER="/apps/Miniconda3-latest-Linux-x86_64.sh"

    # Set LC_ALL=C to avoid possible errors when install packages via pip
    # Reference: https://stackoverflow.com/questions/36394101/pip-install-locale-error-unsupported-locale-setting
    export LC_ALL=C

    # check if conda has already been installed
    conda --version >/dev/null 2>&1
    install_conda=$?

    # download and install miniconda
    if [ $install_conda -ne 0 ]; then
        # default installation directory
        INSTALLATION_DIR=$HOME/storage/miniconda
        echo "-- Downloading Environment manager (conda) --"
        echo "After that, its path will be added into your .bashrc file."
        echo "conda's path: $INSTALLATION_DIR"

        if [ -d $INSTALLATION_DIR ]; then
            echo -e "\n$INSTALLATION_DIR directory exists. To continue, this directory will be removed."
            read -p "Do you want to continue? (y/n) " answer
            if echo "$answer" | grep -iq "^y"; then
                rm -r $INSTALLATION_DIR
            else
                exit
            fi
        fi

        if [ ! -f $INSTALLER ]; then
            wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
            bash miniconda.sh -b -p $INSTALLATION_DIR
        else
            bash $INSTALLER -b -p $INSTALLATION_DIR
        fi
        echo "Existing ~/.bashrc willl be backed up as ~/.bashrc.python_bak"
        rsync -az ~/.bashrc{,.python_bak}
        echo "export PATH=$INSTALLATION_DIR/bin:"'$PATH' >>~/.bashrc
        source ~/.bashrc
    else
        echo "conda exists at $(which conda). Skip installing conda..."
    fi
    # note that we need to update conda to 4.10+
    echo "Checking installed conda version...."
    conda --version
    read -p "Conda version is required to be 4.10+. Do you want to update conda? (y/n) " answer
    if echo "$answer" | grep -iq "^y"; then
        conda update conda -c conda-forge
    fi

    # the default environment is python 3.6
    echo "-- Setting up Python packages --"
    echo "The default Python environment to be set up is Python 3.6."
    # check whether user has installed CBIG_py3 previously
    cbig_py3_binfile=$HOME/storage/miniconda/envs/CBIG_py3/bin/python3
    cbig_py3_bak_dir=$HOME/storage/miniconda/envs/CBIG_py3_bak
    echo "-- Installing Python 3.6 --"
    if [ -f $cbig_py3_binfile ]; then
        echo -e "CBIG_py3 already exists. To continue, current CBIG_py3 would be backed up as CBIG_py3_bak."
        read -p "Do you want to continue? (y/n) " answer
        if echo "$answer" | grep -iq "^y"; then
            echo "-- Backing up CBIG_py3 as CBIG_py3_bak --"
            # check whether cbig_py3_bak_dir already exists
            if [ -d $cbig_py3_bak_dir ]; then
                echo -e "\n$cbig_py3_bak_dir already exists. To continue, this directory will be removed."
                read -p "Do you want to continue? (y/n) " answer
                if echo "$answer" | grep -iq "^y"; then
                    rm -r $cbig_py3_bak_dir
                else
                    exit
                fi
            fi
            conda create --name CBIG_py3_bak --clone CBIG_py3 --yes
            echo "-- Removing current CBIG_py3 --"
            conda remove --name CBIG_py3 --all --yes
        else
            exit
        fi
    fi
    # create new CBIG_py3
    echo "-- Creating CBIG_py3 environment --"
    conda env create --file $WORKDIR/CBIG_python3_env.yml
    source activate CBIG_py3

    conda_failed=$?

    if [ $conda_failed -eq 0 ]; then
        echo "Success!"
        if [ $install_conda -ne 0 ]; then
            echo "Please log out and log in again to complete the installation"
        fi
    else
        echo "Failed to install"
        if [ $install_conda -ne 0 ]; then
            rm -r $INSTALLATION_DIR
            rm ~/.bashrc
            mv ~/.bashrc{.python_bak,}
            echo "Reverted back to original state"
        fi
    fi
else
    echo "Aborting installation..."
fi