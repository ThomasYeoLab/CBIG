#!/bin/bash
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

bold_file=$1
outdir=$2
tmp_dir=$3
low_f=$4
high_f=$5

root_dir=`dirname "$(readlink -f "$0")"`
if [ -z "$bold_file" ] || [ -z "$outdir" ] || [ -z "$outdir" ]; then
    echo "bold_file, outdir, or tmp_dir cannot be empty" 1>&2
    exit 1
fi

if [ -z "$low_f" ];then
    echo "running fsl_motion_outliers to compute FDRMS"
    fsl_motion_outliers -i $bold_file -o $outdir/${bold_file}_motion_outliers_confound_FDRMS \
                        -s $outdir/${bold_file}_motion_outliers_FDRMS -p $outdir/${bold_file}_motion_outliers_FDRMS \
                        -t $tmp_dir --fdrms
else
    echo "running respiratory pseudomotion filtering"
    mkdir -p $tmp_dir
    refnum=`$FSLDIR/bin/fslval $bold_file dim4`;
    refnum=`echo $refnum / 2 | bc`;
    file_base=`basename $bold_file`
    $FSLDIR/bin/mcflirt -in $bold_file -out $tmp_dir/$file_base -mats -plots -refvol $refnum -rmsrel -rmsabs
    motion_file=$tmp_dir/${file_base}.par
    motion_regressors=${bold_file}_mc.par
    save_FD=$outdir/${file_base}_motion_outliers_FDRMS
    matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; CBIG_preproc_filter_out_respiratory_pseudomotion( \
    '$bold_file', '$motion_file', '$motion_regressors', '$save_FD', '$low_f', '$high_f'); exit; "
    rm -rf $tmp_dir/
fi

exit 0
