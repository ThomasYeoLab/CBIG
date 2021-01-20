#!/bin/sh
# This function replicates the variance component model results in the GSP dataset shown in Li et al., 2019
# Only two behavioral measures are included: Shipley_Vocab_Raw and Matrix_WAIS
#
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

########################
# setup for CIRC cluster
########################
curr_dir=$(pwd)
work_dir=$HOME/cluster/

echo $curr_dir
echo $work_dir

if [ ! -d $work_dir ]; then
	mkdir -p $work_dir
fi

cd $work_dir

########################
# common input variables
########################
project_dir="$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR"
replication_dir="$project_dir/unit_tests/intelligence_score"

test_dir=$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/VarianceComponentModel/GSP
subject_list="$test_dir/lists/subject_list_862.txt"
FD_file="$test_dir/lists/FD_regressor_862.txt"
DVARS_file="$test_dir/lists/DV_regressor_862.txt"
data_csv="$test_dir/lists/GSP_phenotypes_shuffled.csv"
d=431
num_samples=5
rmsub_prefix="subjects862"

top_outdir=$1

for pipeline in GSR Baseline ; do
	RSFC_file=$test_dir/RSFC_862_Fisher_${pipeline}.mat
	outdir=$top_outdir/$pipeline
	
	##########################
	# 2 behavioral measures
	##########################
	cog_list="$replication_dir/scripts/GSP_lists/intelligence_score.txt"
	covariate_list="$replication_dir/scripts/GSP_lists/covariates.txt"
	ystem=2intelligence
	
	log_file="${top_outdir}/CBIG_LiGSR_LME_unittest_intelligence_score_GSP_${pipeline}.log"
	cmd="$project_dir/VarianceComponentModel/scripts/CBIG_LiGSR_LME_workflowGSP.sh -RSFC_file $RSFC_file -trait_list "
	cmd="$cmd $cog_list -covariate_list $covariate_list -FD_file $FD_file -DVARS_file $DVARS_file -subject_list "
	cmd="$cmd $subject_list -outdir $outdir -ystem $ystem -d $d -num_samples $num_samples -rmsub_prefix $rmsub_prefix "
	cmd="$cmd -data_csv $data_csv"
    cmd="$cmd | tee -a ${log_file}"

    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 1:00:00 -mem 8G \
    -name "LiGSRUT_ME"

	sleep 3s
done


