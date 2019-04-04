#!/bin/sh
# This function uses kernel ridge regression to predict fluid intelligence score (PMAT24_A_CR) in the HCP dataset.
#
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

########################
# setup for CIRC cluster
########################
curr_dir=$(pwd)
username=$(whoami)
work_dir=/data/users/$username/cluster/

echo $curr_dir
echo $username
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

test_dir=/mnt/eql/yeo1/CBIG_private_data/unit_tests/stable_projects/preprocessing/Li2019_GSR/intelligence_score/\
KernelRidgeRegression/HCP
subject_list="$test_dir/lists/subject_list_953.txt"
FD_file="$test_dir/lists/FD_regressor_953.txt"
DVARS_file="$test_dir/lists/DV_regressor_953.txt"
#RSFC_file="$test_dir/cort+subcort_new_S1200_953_Fisher.mat"

top_outdir=$1
#top_outdir=$test_dir/ref_output

for pipeline in GSR Baseline; do
	RSFC_file="$test_dir/cort+subcort_new_S1200_953_Fisher_${pipeline}.mat"
	outdir=$top_outdir/$pipeline
	
	##########################
	# call the GSP wrapper
	##########################
	cog_list="$replication_dir/scripts/HCP_lists/PMAT24.txt"
	covariate_list="$replication_dir/scripts/HCP_lists/covariates.txt"
	outstem=PMAT24
	
	for seed in $(seq 1 1 3); do
		cmd="$project_dir/KernelRidgeRegression/HCP/scripts/CBIG_LiGSR_KRR_workflowHCP.sh -subject_list $subject_list "
		cmd="$cmd -RSFC_file $RSFC_file -y_list $cog_list -covariate_list $covariate_list -FD_file $FD_file -DVARS_file "
		cmd="$cmd $DVARS_file -outdir $outdir -outstem $outstem -seed $seed -num_test_folds 5 -num_inner_folds 5"
		
		echo $cmd | qsub -V -q circ-spool -l walltime=2:00:00,mem=6GB -m ae -N CBIG_LiGSR_KRR_unittest_PMAT_HCP
		
		if [ ! -f $outdir/covariates_${outstem}.mat ] || [ ! -f $outdir/y_${outstem}.mat ]; then
			# wait for the files shared across random splits to be saved
			sleep 2m   
		else
			sleep 3s
		fi
	done
done


