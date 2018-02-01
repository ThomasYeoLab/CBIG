#!/bin/csh
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set orig_data_dir = "$CBIG_CODE_DIR/data/example_data"
set output_dir = $1
set sub_list = "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/examples/scripts/example_sub_list.txt"
set subjects = `cat $sub_list`
set code_dir = "${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering"
set surf_stem = "_rest_skip4_stc_mc_resid_cen_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5"
set outlier_stem = "_FDRMS0.2_DVARS50_motion_outliers"

## Create folder structure within output_dir, and make soft links of input files to orig_data_dir
foreach s ($subjects)
	set s_id = `echo $s | cut -d '_' -f 1`
	set sess_id = `echo $s | cut -d '_' -f 2`
	mkdir -p $output_dir/subjects/$s
	ln -s $orig_data_dir/$s_id/$s/qc $output_dir/subjects/$s/
	ln -s $orig_data_dir/$s_id/$s/surf $output_dir/subjects/$s/
	ln -s $orig_data_dir/$s_id/$s/logs $output_dir/subjects/$s/
end


## Call wrapper function
mkdir -p $output_dir/clustering
${code_dir}/CBIG_Yeo2011_general_cluster_fcMRI_surf2surf_profiles.csh -sd ${output_dir}/subjects -sub_ls ${sub_list} -surf_stem ${surf_stem} -n 17 -out_dir ${output_dir}/clustering -cluster_out ${output_dir}/clustering/HNU_example_clusters017_scrub -tries 5 -outlier_stem ${outlier_stem}
