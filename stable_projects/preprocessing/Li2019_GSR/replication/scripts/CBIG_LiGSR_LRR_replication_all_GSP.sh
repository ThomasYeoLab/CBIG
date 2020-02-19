#!/bin/sh
# This function replicate the linear ridge regression results in the GSP dataset shown in Li et al., 2019
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
replication_dir="$project_dir/replication"

test_dir=$CBIG_REPDATA_DIR/stable_projects/preprocessing/Li2019_GSR/LinearRidgeRegression/GSP
subject_list="$test_dir/lists/subject_list_862.txt"
FD_file="$test_dir/lists/FD_regressor_862.txt"
DVARS_file="$test_dir/lists/DV_regressor_862.txt"

gpso_dir=$CBIG_CODE_DIR/external_packages/matlab/non_default_packages/Gaussian_Process
top_outdir=$1

for pipeline in GSR Baseline ; do
	RSFC_file=$test_dir/RSFC_862_Fisher_${pipeline}.mat
	outdir=$top_outdir/$pipeline
	
	##########################
	# 23 behavioral measures
	##########################
	cog_list="$replication_dir/scripts/GSP_lists/23behaviors.txt"
	covariate_list="$replication_dir/scripts/GSP_lists/covariates_23behaviors.txt"
	y_names=$(cat $cog_list)
	
	for y_name in $y_names; do
		for seed in $(seq 1 1 20); do
			cmd="$project_dir/LinearRidgeRegression/GSP/scripts/CBIG_LiGSR_LRR_workflowGSP.sh -subject_list "
			cmd="$cmd $subject_list -RSFC_file $RSFC_file -y_name $y_name -covariate_list $covariate_list -FD_file "
			cmd="$cmd $FD_file -DVARS_file $DVARS_file -outdir $outdir -gpso_dir $gpso_dir -seed $seed "
			
			echo $cmd | $CBIG_SCHEDULER_DIR/qsub -V -q circ-spool -l walltime=21:00:00,mem=7GB -m ae \
			  -N CBIG_LiGSR_LRR_replication_all_GSP
			
			if [ ! -f $outdir/covariates/covariates.mat ] || [ ! -f $outdir/y/y_${y_name}.mat ]; then
				# wait for the files shared across random splits to be saved
				sleep 3m   
			else
				sleep 3s
			fi
		done
	done
	
	
	##########################
	# predict age
	##########################
	age_list="$replication_dir/scripts/GSP_lists/Age_header.txt"
	covariate_list="$replication_dir/scripts/GSP_lists/covariates_age.txt"
	y_name=$(cat $age_list)
	
	for seed in $(seq 1 1 20); do
		cmd="$project_dir/LinearRidgeRegression/GSP/scripts/CBIG_LiGSR_LRR_workflowGSP.sh -subject_list $subject_list "
		cmd="$cmd -RSFC_file $RSFC_file -y_name $y_name -covariate_list $covariate_list -FD_file $FD_file -DVARS_file "
		cmd="$cmd $DVARS_file -outdir $outdir -gpso_dir $gpso_dir -seed $seed "
		
		echo $cmd | $CBIG_SCHEDULER_DIR/qsub -V -q circ-spool -l walltime=21:00:00,mem=7GB -m ae \
		  -N CBIG_LiGSR_LRR_replication_all_GSP
		
		if [ ! -f $outdir/covariates/covariates_${y_name}.mat ] || [ ! -f $outdir/y/y_${y_name}.mat ]; then
			sleep 3m
		else
			sleep 3s
		fi
	done
	
	
done


