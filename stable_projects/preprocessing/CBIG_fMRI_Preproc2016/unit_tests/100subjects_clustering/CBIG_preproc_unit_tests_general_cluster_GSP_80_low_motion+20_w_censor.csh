#!/bin/csh -f
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set your_out_dir = $1  # where you want to put the clustering outputs
set data_dir = $2    # where the 100 subjects preprocessed by you are stored
set subject_list = "$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests"
set subject_list = "$subject_list/100subjects_clustering/GSP_80_low_motion+20_w_censor.txt"

# create fake dir for each subject in your data dir, and make symbolic links to original data dir
set your_subject_dir = "${your_out_dir}/preproc_out"
foreach s (`cat $subject_list`)
    mkdir -p ${your_subject_dir}/$s
    ln -s ${data_dir}/$s/surf ${your_subject_dir}/$s/
    ln -s ${data_dir}/$s/logs ${your_subject_dir}/$s/
    ln -s ${data_dir}/$s/qc ${your_subject_dir}/$s/
    ln -s ${data_dir}/$s/bold ${your_subject_dir}/$s/
end

set out_dir = "${your_out_dir}/clustering"
set code_dir = "${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering"
set num_clusters = 17
set formated_cluster = `echo $num_clusters | awk '{printf ("%03d", $1)}'`
set scrub_flag = 1   # 1 for scrub, 0 for not scrub

set surf_stem = "_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5"
set outlier_stem = "_FDRMS0.2_DVARS50_motion_outliers"

set curr_dir = `pwd`
set work_dir = $HOME/cluster/ 

echo $curr_dir
echo $work_dir

if (! -e $work_dir) then
    mkdir -p $work_dir
endif

cd $work_dir

set log_file = "${out_dir}/clust_100sub_ut.log"

if( $scrub_flag == 1 ) then
    set cmd = "${code_dir}/CBIG_Yeo2011_general_cluster_fcMRI_surf2surf_profiles.csh -sd"
    set cmd = "$cmd ${your_subject_dir} -sub_ls ${subject_list} -surf_stem ${surf_stem} -n ${num_clusters} -out_dir"
    set cmd = "$cmd ${out_dir} -cluster_out ${out_dir}/GSP_80_low_mt_20_w_censor_clusters${formated_cluster}_scrub"
    set cmd = "$cmd -tries 1000 -outlier_stem ${outlier_stem}"
else
    set cmd = "${code_dir}/CBIG_Yeo2011_general_cluster_fcMRI_surf2surf_profiles.csh -sd" 
    set cmd = "$cmd ${your_subject_dir} -sub_ls ${subject_list} -surf_stem ${surf_stem} -n ${num_clusters} -out_dir"
    set cmd = "$cmd ${out_dir} -cluster_out \
          ${out_dir}/GSP_80_low_mt_20_w_censor_clusters${formated_cluster}_noscrub"
    set cmd = "$cmd -tries 1000 "
endif

set cmd = "$cmd | tee -a ${log_file}"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 20:00:00 -mem 2G -name "clust_100sub_ut" 



