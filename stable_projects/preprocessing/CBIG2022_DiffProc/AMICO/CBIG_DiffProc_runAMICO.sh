#!/bin/sh
#####
# Example: 
#    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_diffusion_processing2022/AMICO/CBIG_DiffProc_runAMICO.sh \
#        --subj_list /path/to/txtfile --dwi_dir /path/to/dwi_images \
#        --output_dir /path/to/output --py_env name_of_AMICO_environment \
#        --mask_output_dir /path/to/stored/b0_brainmask
#
# This function fits the NODDI model for a given list of subjects using the AMICO pipeline.
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

# set up directories and lists
scriptdir=$( dirname "${BASH_SOURCE[0]}" )

while [[ $# -gt 0 ]]; do
key="$1"

case $key in

    -s|--subj_list)
    subj_list="$2"
        shift; shift;;

    -d|--dwi_dir)
    dwi_dir="$2/"
        shift; shift;;

    -o|--output_dir)
    output_dir="$2/"
        shift; shift;;

    -p|--py_env)
    py_env="$2"
        shift; shift;;

    -m|--mask_output_dir)
    mask_output_dir="$2/"
        shift; shift;;

    *)    # unknown option
    echo "Unknown option: $1"
    shift;;
esac
done

# check args
if [ -z "$subj_list" ] || [ -z "$dwi_dir" ] || [ -z "$output_dir" ] || [ -z "$py_env" ]; then
    echo "Missing compulsory variable!"
    exit
fi

if [ -z "$mask_output_dir" ]; then
    mask_output_dir="NIL"
fi

# print selected options
echo "[subj_list] = $subj_list"
echo "[dwi dir] = $dwi_dir"
echo "[output dir] = $output_dir"
echo "[py env] = $py_env"
echo "[mask_output_dir] = $mask_output_dir"

# create directories
if [ ! -d  $output_dir/output ]; then mkdir -p $output_dir/output; fi
log_dir=$output_dir/logs 
if [ ! -d  $log_dir ]; then mkdir -p $log_dir; fi

# run for all subjects
cat $subj_list| while read subj; do
    echo " ------ Processing $subj ------ "
    if [ ! -e $output_dir/output/AMICO/NODDI/$subj/FIT_OD.nii.gz ]; then # skip if output file exists

        # generate b0 mask
        if [[ $mask_output_dir == 'NIL' ]]; then
            mask_output_dir=$output_dir/b0_mask
            if [ ! -d  $mask_output_dir/$subj ]; then mkdir -p $mask_output_dir/$subj; fi
            echo "Generating brain mask from b=0 image..."
            if [ ! -e  $mask_output_dir/$subj/${subj}_bet_b0_mask.nii.gz ]; then
                $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG2022_DiffProc/utilities/CBIG_DiffProc_generate_b0mask.sh \
                    ${dwi_dir}/${subj} $subj $mask_output_dir
            else
                echo "Brain mask already exists, skipping..."
            fi
        fi
        mask=$mask_output_dir/${subj}/${subj}_bet_b0_mask.nii.gz

        # run AMICO
        cmd="source activate $py_env; python $scriptdir/CBIG_DiffProc_runAMICO.py \
            $output_dir $subj $dwi_dir $mask > $log_dir/AMICO_${subj}_log.txt"
        ssh headnode "$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd '$cmd' -walltime '5:00:00' -mem '15G' \
            -name 'AMICO' -joberr '$log_dir' -jobout '$log_dir' " < /dev/null

        # wait for kernel to generate before running next job
        if [ ! -d  $output_dir/AMICO/kernels ]; then
            echo "Waiting for kernel to generate..."
            sleep 5m
        fi
    else
        echo "$subj exists, skipping..."
    fi
done
