#!/bin/bash
#####
# Example: 
#    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG2022_DiffProc/ \
#        MRtrix/CBIG_DiffProc_batch_tractography.sh \
#        --subj_list /path/to/txtfile --dwi_dir /path/to/dwi_images \
#        --output_dir /path/to/output --py_env name_of_AMICO_environment \
#        --mask_output_dir /path/to/stored/b0_brainmask
#
# This function fits a tractogram and generates a structural connectivity matrix for each subject
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

###############
# set up environment
###############

# set up directories and lists
scriptdir=$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG2022_DiffProc/MRtrix/

while [[ $# -gt 0 ]]; do
key="$1"

case $key in

    -s|--subj_list)
    subj_list="$2"
        shift; shift;;

    -i|--input_dir)
    input_dir="$2"
        shift; shift;;

    -o|--output_dir)
    output_dir="$2"
        shift; shift;;

    -p|--py_env)
    py_env="$2"
        shift; shift;;

    -m|--mask_output_dir)
    mask_output_dir="$2"
        shift; shift;;

    -a|--algo)
    algo="$2"
        shift; shift;;

    -t|--tract_streamlines)
    tract_streamlines="$2"
        shift; shift;;

    *)    # unknown option
    echo "Unknown option: $1"
    shift;;
esac
done

# check args
if [ -z "$subj_list" ] || [ -z "$input_dir" ] || [ -z "$output_dir" ] || [ -z "$py_env" ]; then
    echo "Missing compulsory variable!"
    exit
fi

if [ -z "$mask_output_dir" ]; then
    mask_output_dir="NIL"
fi

if [ -z "$algo" ]; then
    algo='iFOD2'
fi

if [ -z "$tract_streamlines" ]; then
    tract_streamlines="5M"
fi

# print selected options
echo "[subj_list] = $subj_list"
echo "[input dir] = $input_dir"
echo "[output dir] = $output_dir"
echo "[py env] = $py_env"
echo "[mask_output_dir] = $mask_output_dir"
echo "[algo] = $algo"
echo "[tract_streamlines] = $tract_streamlines"

logdir=$output_dir/$algo/logs
if [ ! -d $logdir ]; then mkdir -p $logdir; fi

###############
# 1. generate required files
###############
# set up list of parcellations
parcellation_arr="Schaefer2018_100Parcels_17Networks,Schaefer2018_400Parcels_17Networks"
parcels_num_arr="100,400"

cat $subj_list | while read subject; do
    # if subject name needs to be modified (e.g. remove prefix), it should be done here
    # for example, delete underscore : sub=$( echo $subject | tr -d '_')
    sub=$( echo $subject )
    echo "[PROGRESS]: Generating parcellations and tractogram for: $sub"
    $scriptdir/CBIG_DiffProc_preprocess_tractography.sh $sub $algo $scriptdir $input_dir $output_dir $py_env \
        $mask_output_dir $tract_streamlines $parcellation_arr $parcels_num_arr
done

###############
# 1.5 wait for job completion
###############
echo "[PROGRESS]: Waiting for all jobs to complete before moving on..."
start=$SECONDS
user_id=$( whoami )
parc_jobs=$( ssh headnode "qselect -u $user_id -N gen_parcellation | wc -l" )
trac_jobs=$( ssh headnode "qselect -u $user_id -N gen_tractogram | wc -l" )
while [[ $parc_jobs -gt 0 ]] || [[ $trac_jobs -gt 0 ]]; do
    duration_min=$((( $SECONDS - $start ) / 60))
    echo "    ${duration_min}m: $parc_jobs parcellation jobs & $trac_jobs tractogram jobs still running"
    sleep 3m
    parc_jobs=$( ssh headnode "qselect -u $user_id -N gen_parcellation | wc -l" )
    trac_jobs=$( ssh headnode "qselect -u $user_id -N gen_tractogram | wc -l" )
done

###############
# 2. generate connectome
###############
cat $subj_list | while read subject; do
    # if subject name needs to be modified (e.g. remove prefix), it should be done here
    # for example, delete underscore : sub=$( echo $subject | tr -d '_')
    sub=$( echo $subject )
    echo "[PROGRESS]: Generating connectome for subj ID: $sub"
    cmd="$scriptdir/CBIG_DiffProc_tractography_3_gen_connectome.sh $sub $algo $input_dir $output_dir \
        $tract_streamlines $parcellation_arr $parcels_num_arr > $logdir/${sub}_SC_log.txt 2>&1"
    ssh headnode "source activate $py_env; \
    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd '$cmd' -walltime 02:30:00 \
    -name 'gen_connectome' -mem 6GB -joberr '$logdir' -jobout '$logdir'" < /dev/null
done

###############
# 2.5 wait for job completion
###############
echo "[PROGRESS]: Waiting for all jobs to complete before moving on..."
start=$SECONDS
user_id=$( whoami )
conn_jobs=$( ssh headnode "qselect -u $user_id -N gen_connectome | wc -l" )
while [[ $conn_jobs -gt 0 ]]; do
    duration_min=$((( $SECONDS - $start ) / 60))
    echo "    ${duration_min}m: $conn_jobs connectome jobs still running"
    sleep 3m
    conn_jobs=$( ssh headnode "qselect -u $user_id -N gen_connectome | wc -l" )
done

###############
# 3. check logs and replace missing parcels with NaN
###############
cat $subj_list | while read subject; do
    # if subject name needs to be modified (e.g. remove prefix), it should be done here
    # for example, delete underscore : sub=$( echo $subject | tr -d '_')
    sub=$( echo $subject )
    echo "[PROGRESS]: Checking connectomes and filling missing nodes with NaN for subj ID: $sub"
    matlab -nodesktop -nosplash -nodisplay -r " try addpath('$scriptdir'); \
        CBIG_DiffProc_mrtrix_fill_na('$logdir/${sub}_SC_log.txt'); catch ME; \
        display(ME.message); end; exit; " >> $logdir/${sub}_SC_log.txt 2>&1
done
