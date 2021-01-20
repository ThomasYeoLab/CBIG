#!/bin/sh
# This function replicate the linear ridge regression results in the GSP dataset shown in Li et al., 2019
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

test_dir=$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/LinearRidgeRegression/GSP
subject_list="$test_dir/lists/subject_list_862.txt"
FD_file="$test_dir/lists/FD_regressor_862.txt"
DVARS_file="$test_dir/lists/DV_regressor_862.txt"
data_csv="$test_dir/lists/GSP_phenotypes_shuffled.csv"

gpso_dir=$1
top_outdir=$2
evaluations=5
tree=2

for pipeline in GSR Baseline ; do
	RSFC_file=$test_dir/RSFC_862_Fisher_${pipeline}.mat
	outdir=$top_outdir/$pipeline
	
	##########################
	# call the GSP wrapper
	##########################
	cog_list="$replication_dir/scripts/GSP_lists/intelligence_score.txt"
	covariate_list="$replication_dir/scripts/GSP_lists/covariates.txt"
	y_names=$(cat $cog_list)
	
	for y_name in $y_names; do
		for seed in $(seq 1 1 3); do
		    log_file="${top_outdir}/CBIG_LiGSR_LRR_unittest_intelligence_score_GSP_${y_name}_${seed}.log"
			cmd="$project_dir/LinearRidgeRegression/GSP/scripts/CBIG_LiGSR_LRR_workflowGSP.sh -subject_list "
			cmd="$cmd $subject_list -RSFC_file $RSFC_file -y_name $y_name -covariate_list $covariate_list -FD_file "
			cmd="$cmd $FD_file -DVARS_file $DVARS_file -outdir $outdir -gpso_dir $gpso_dir -seed $seed -num_test_folds "
			cmd="$cmd 3 -num_inner_folds 3 -eval $evaluations -tree $tree -data_csv $data_csv"
            cmd="$cmd | tee -a ${log_file}"

            $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 1:00:00 -mem 7G \
            -name "LiGSRUT_LR"
			
			if [ ! -f $outdir/covariates/covariates.mat ] || [ ! -f $outdir/y/y_${y_name}.mat ]; then
				# wait for the files shared across random splits to be saved
				sleep 1m   
			else
				sleep 3s
			fi
		done
	done
	
done


