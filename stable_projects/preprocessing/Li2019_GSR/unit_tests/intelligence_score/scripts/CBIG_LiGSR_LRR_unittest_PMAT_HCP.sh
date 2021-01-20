#!/bin/sh
# This function uses linear ridge regression to predict fluid intelligence score (PMAT24_A_CR) in the HCP dataset.
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

test_dir=$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/LinearRidgeRegression/HCP
subject_list="$test_dir/lists/fake_subject_list_953.txt"
FD_file="$test_dir/lists/FD_regressor_953.txt"
DVARS_file="$test_dir/lists/DV_regressor_953.txt"
restricted_csv="$test_dir/lists/fake_family.csv"
unrestricted_csv="$test_dir/lists/fake_behavior.csv"

gpso_dir=$1
top_outdir=$2
evaluations=5
tree=2

for pipeline in GSR Baseline; do
	RSFC_file="$test_dir/cort+subcort_new_S1200_953_Fisher_${pipeline}.mat"
	outdir=$top_outdir/$pipeline
	
	##########################
	# call the GSP wrapper
	##########################
	cog_list="$replication_dir/scripts/HCP_lists/PMAT24.txt"
	covariate_list="$replication_dir/scripts/HCP_lists/covariates.txt"
	y_name=$(cat $cog_list)
	
	for seed in $(seq 1 1 3); do
		log_file="${top_outdir}/CBIG_LiGSR_LRR_unittest_PMAT_HCP_${seed}.log"
		cmd="$project_dir/LinearRidgeRegression/HCP/scripts/CBIG_LiGSR_LRR_workflowHCP.sh -subject_list $subject_list "
		cmd="$cmd -RSFC_file $RSFC_file -y_name $y_name -covariate_list $covariate_list -FD_file $FD_file -DVARS_file "
		cmd="$cmd $DVARS_file -outdir $outdir -gpso_dir $gpso_dir -seed $seed -num_test_folds 3 -num_inner_folds 3 "
		cmd="$cmd -restricted_csv $restricted_csv -unrestricted_csv $unrestricted_csv -eval $evaluations -tree $tree"
        cmd="$cmd | tee -a ${log_file}"

        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 2:00:00 -mem 6G \
        -name "LiGSRUT_LR"
		
		if [ ! -f $outdir/covariates/covariates.mat ] || [ ! -f $outdir/y/y_${y_name}.mat ]; then
			# wait for the files shared across random splits to be saved
			sleep 2m   
		else
			sleep 3s
		fi
	done
done


