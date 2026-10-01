#!/bin/bash
#
# This wrapper script computes Framewise Displacement (FD) or RMS-based FD (FDRMS)
# for a given BOLD fMRI run.
#
# If no frequencies are specified (e.g., no `low_f`), the script directly computes FD/FDRMS using
# FSL's fsl_motion_outliers.
#
# If pseudorespiratory low-pass (`low_f`) and/or high-pass (`high_f`) frequencies are provided,
# the script performs respiratory pseudomotion filtering using MATLAB
# (CBIG_preproc_filter_out_respiratory_pseudomotion), then computes the
# filtered FD metric.
# If only `low_f` provided (and `high_f` not provided), low-pass filtering is performed.
# But if both `low_f` and `high_f` are provided, bandstop filtering is performed.
#
# Example:
# CBIG_preproc_compute_FD_metric_wrapper.sh -i "sub-01_bld001_rest_skip8_stc" \
#     -o "base_outdir/sub-01/bold/mc" -t "base_outdir/sub-01/bold/mc/tmp_outliers/001" \
#     -f "FD" -l "0.2" -u "0.33" \
#     -m "$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/matlab_runtime/utilities"
#
# Written by Jianzhong Chen, Trevor Tan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Initialize variables
bold_file=""
outdir=""
tmp_dir=""
FD_metric=""
low_f=""
high_f=""
matlab_runtime_util=""
root_dir=`dirname "$(readlink -f "$0")"`

# Parse options
while getopts "i:o:t:f:l:u:m:" arg; do
    case $arg in
        i) bold_file=$OPTARG;;
        o) outdir=$OPTARG;;
        t) tmp_dir=$OPTARG;;
        f) FD_metric=$OPTARG;;
        l) low_f=$OPTARG;;
        u) high_f=$OPTARG;;
        m) matlab_runtime_util=$OPTARG;;
        ?)
        echo "Unknown argument"
        exit 1;;
    esac
done

# Check arguments
if [ -z "$bold_file" ] || [ -z "$outdir" ] || [ -z "$tmp_dir" ] || [ -z "$FD_metric" ]; then
    echo "bold_file, outdir, or tmp_dir cannot be empty" 1>&2
    exit 1
fi

if [[ "$FD_metric" != "FD" && "$FD_metric" != "FDRMS" ]]; then
    echo "ERROR: FD_metric can only be FDRMS or FD." 1>&2
    exit 1
fi

if [ ! -d "$tmp_dir" ]; then
    mkdir -p "$tmp_dir"
fi

# If no frequencies provided, no respiratory pseudomotion filtering performed, run FSL function.
# But if provided, run respiratory pseudomotion filtering.
if [ -z "$low_f" ];then
    if [ "$FD_metric" == "FDRMS" ]; then
        echo "Running fsl_motion_outliers to compute FDRMS"
        cmd="fsl_motion_outliers -i $bold_file -o $outdir/${bold_file}_motion_outliers_confound_FDRMS"
        cmd="$cmd -s $outdir/${bold_file}_motion_outliers_FDRMS -p $outdir/${bold_file}_motion_outliers_FDRMS"
        cmd="$cmd -t $tmp_dir --fdrms"
    else
        echo "Running fsl_motion_outliers to compute FD"
        cmd="fsl_motion_outliers -i $bold_file -o $outdir/${bold_file}_motion_outliers_confound_FD"
        cmd="$cmd -s $outdir/${bold_file}_motion_outliers_FD -p $outdir/${bold_file}_motion_outliers_FD"
        cmd="$cmd -t $tmp_dir --fd"
    fi
    echo "$cmd"
    eval "$cmd"
else
    echo "Running respiratory pseudomotion filtering"
    refnum=`$FSLDIR/bin/fslval $bold_file dim4`;
    refnum=`echo $refnum / 2 | bc`;
    file_base=`basename $bold_file`
    $FSLDIR/bin/mcflirt -in $bold_file -out $tmp_dir/$file_base -mats -plots -refvol $refnum -rmsrel -rmsabs
    motion_file=$tmp_dir/${file_base}.par
    motion_regressors=${bold_file}_mc.par
    save_FD=$outdir/${file_base}_motion_outliers_${FD_metric}
    matlab_args_array=("'$bold_file'" "'$motion_file'" "'$motion_regressors'" "'$FD_metric'" "'$save_FD'"
        "'$low_f'" "'$high_f'")

    # If no path to MATLAB Runtime utilities provided, run using MATLAB
    # But if provided, run using MATLAB Runtime
    if [ -z "$matlab_runtime_util" ]; then
        echo "Running MATLAB version for CBIG_preproc_filter_out_respiratory_pseudomotion"
        # Insert commas between arguments
        matlab_args=""
        for element in "${matlab_args_array[@]}"; do
            matlab_args="${matlab_args}${element},"
        done
        matlab_args="${matlab_args%,}" # remove trailing comma
        cmd="matlab -batch \"addpath('${root_dir}');"
        cmd="$cmd CBIG_preproc_filter_out_respiratory_pseudomotion($matlab_args);\""
        echo "$cmd"
        eval "$cmd"
    else
        echo "Running MATLAB Runtime version for CBIG_preproc_filter_out_respiratory_pseudomotion"
        cmd="${matlab_runtime_util}/CBIG_preproc_filter_out_respiratory_pseudomotion ${matlab_args_array[@]}"
        echo "$cmd"
        eval "$cmd"
    fi
    rm -rf $tmp_dir/
fi

exit 0
