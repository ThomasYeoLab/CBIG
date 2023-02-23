#! /bin/bash

# Written by Tianchu Zeng and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

echo "Executing this script will install a new environment called CBIG_py3_aws and install package aws-cli."
read -p "Are you sure? (y/n) " answer

if echo "$answer" | grep -iq "^y"; then
    echo "The default Python environment to be set up is Python 3.6."
    # check whether user has installed CBIG_py3_aws previously
    cbig_py3_aws_binfile=$HOME/storage/miniconda/envs/CBIG_py3_aws/bin/python3
    cbig_py3_aws_bak_dir=$HOME/storage/miniconda/envs/CBIG_py3_aws_bak
    echo "-- Installing Python 3.6 --"
    if [ -f $cbig_py3_aws_binfile ]; then
        echo -e "CBIG_py3_aws already exists. To continue, current CBIG_py3_aws would be backed up as CBIG_py3_aws_bak."
        read -p "Do you want to continue? (y/n) " answer
        if echo "$answer" | grep -iq "^y"; then
            echo "-- Backing up CBIG_py3_aws as CBIG_py3_aws_bak --"
            # check whether cbig_py3_bak_dir already exists
            if [ -d $cbig_py3_aws_bak_dir ]; then
                echo -e "\n$cbig_py3_aws_bak_dir already exists. To continue, this directory will be removed."
                read -p "Do you want to continue? (y/n) " answer
                if echo "$answer" | grep -iq "^y"; then
                    rm -r $cbig_py3_aws_bak_dir
                else
                    exit
                fi
            fi
            conda create --name CBIG_py3_aws_bak --clone CBIG_py3_aws --yes
            echo "-- Removing current CBIG_py3_aws --"
            conda remove --name CBIG_py3_aws --all --yes
        else
            exit
        fi
    fi
    echo "-- Creating CBIG_py3_aws environment --"
    conda create -n CBIG_py3_aws python=3.6
    source activate CBIG_py3_aws

    echo "-- Installing package aws-cli --"
    conda install -c conda-forge awscli=1.15.15
    conda deactivate
else
    exit
fi