#! /bin/bash

# Written by Wei Hou Tan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

echo "Executing this script will install Python packages used by CBIG."
read -p "Are you sure? (y/n) " answer

if echo "$answer" | grep -iq "^y" ; then
    # store current directory
    WORKDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    INSTALLER="/apps/Miniconda3-latest-Linux-x86_64.sh"
    
    # Set LC_ALL=C to avoid possible errors when install packages via pip
    # Reference: https://stackoverflow.com/questions/36394101/pip-install-locale-error-unsupported-locale-setting
    export LC_ALL=C

    # check if conda has already been installed
    conda --version > /dev/null 2>&1
    install_conda=$?
    
    # download and install miniconda
    if [ $install_conda -ne 0 ]; then
        # default installation directory
        INSTALLATION_DIR=$HOME/storage/miniconda
        echo "-- Downloading Python packages manager (conda) --"
        echo "After that, its path will be added into your .bashrc file."
        echo "conda's path: $INSTALLATION_DIR"

        if [ -d $INSTALLATION_DIR ]; then
          echo -e "\n$INSTALLATION_DIR directory exists. To continue, this directory will be removed."
          read -p "Do you want to continue? (y/n) " answer
          if echo "$answer" | grep -iq "^y" ; then
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
        echo "export PATH=$INSTALLATION_DIR/bin:"'$PATH' >> ~/.bashrc
        source ~/.bashrc
    else
        echo "conda exists at `which conda`. Skip installing conda..."
    fi
    
    # the default environment is python 3.5
    echo "-- Setting up Python packages --"
    echo "The default Python environment to be set up is Python 3.5."
    read -p "Do you want to continue? (y/n) " answer
    if echo "$answer" | grep -iq "^y" ; then
      echo "-- Installing Python 3.5 --"
      # install python 3 packages into CBIG_py3 environment
      conda create --name CBIG_py3 --file $WORKDIR/CBIG_python3_conda_packages.txt --yes
      source activate CBIG_py3
    else
      read -p "Do you want to set up Python 2.7 instead? (y/n) " answer
      if echo "$answer" | grep -iq "^y" ; then
        echo "-- Installing Python 2.7 --"
        # install python 2 packages into CBIG_py2 environment
        conda create --name CBIG_py2 --file $WORKDIR/CBIG_python2_conda_packages.txt --yes
        source activate CBIG_py2
      else
        exit
      fi
    fi 

    conda_failed=$?
    
    # install common python packages via pip
    # only execute this step if the previous step is successful
    if [ $conda_failed -eq 0 ]; then
        pip install -r $WORKDIR/CBIG_pip_packages.txt
    fi
    pip_failed=$?

    if [[ $conda_failed -eq 0 && $pip_failed -eq 0 ]]; then
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
