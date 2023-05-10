#!/bin/sh
#
# This function check whether the input file is plumb and RPI orientation.
# If not, this function will deoblique the input file and change orientation to RPI.
# This check should always be applied to T1 before recon-all and fMRI before preprocessing.
#
# Example: CBIG_preproc_deoblique.sh -i input.nii.gz
# Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

while getopts "i:o:n" arg
do
    case $arg in
        i) input=$OPTARG;;
        o) output=$OPTARG;;
        n) no_orient="true";;
        ?)
        echo "Unknown argument"
        exit 1;;
    esac
done

# Check whether the input file is oblique
obl=`3dinfo -is_oblique ${input}`

# Check the input file oreintation
ori=`3dinfo -orient ${input}`
if [[ $obl == 1 || $ori != RPI ]]; then

    if [ -n "${output}" ]; then
        3dcopy -overwrite ${input} ${output}
    else
        output=${input}
    fi    
    if [[ $obl == 1 ]]; then
        echo "Input data is oblique. Deoblique to plumb..."
        3drefit -deoblique ${output}
    fi
    if [[ $ori != RPI ]]; then
        if [ ! -n "${no_orient}" ]; then
            echo "Orientation of input data is not RPI. Change to RPI..."
            3dresample -overwrite -orient RPI -inset ${output} -prefix ${output}
        else
            echo "Orientation of input data is not RPI."
        fi
    fi
else
    echo "Input data is already plumb and RPI. No changes required."
    if [ -n "${output}" ]; then
        rsync ${input} ${output}
    fi
fi

