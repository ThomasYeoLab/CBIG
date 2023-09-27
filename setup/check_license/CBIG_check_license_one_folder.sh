#! /bin/sh 

# input should be a directory
path=$1


# for loop all the files and check the license
for file in $path/*
do
    echo $file
    sh $CBIG_CODE_DIR/setup/check_license/CBIG_check_license_matlab_file.sh $file
    clear
done
exit 0;
