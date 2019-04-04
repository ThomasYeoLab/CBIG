#!/bin/sh
# This function replicate the kernel ridge regression results in the HCP dataset shown in Li et al., 2019
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
replication_dir="$project_dir/replication"

test_dir=/mnt/eql/yeo1/CBIG_private_data/replication/stable_projects/preprocessing/Li2019_GSR/KernelRidgeRegression/HCP
subject_list="$test_dir/lists/subject_list_953.txt"
FD_file="$test_dir/lists/FD_regressor_953.txt"
DVARS_file="$test_dir/lists/DV_regressor_953.txt"
#RSFC_file="$test_dir/cort+subcort_new_S1200_953_Fisher.mat"

top_outdir=$1
#top_outdir=$test_dir/ref_output

for pipeline in GSR Baseline ; do
	RSFC_file=$test_dir/cort+subcort_new_S1200_953_Fisher_${pipeline}.mat
	outdir=$top_outdir/$pipeline
	
	##########################
	# 13 cognitive measures
	##########################
	cog_list="$replication_dir/scripts/HCP_lists/Cognitive_unrestricted.txt"
	covariate_list="$replication_dir/scripts/HCP_lists/covariates_58behaviors.txt"
	outstem=13cognitive
	
	for seed in $(seq 1 1 20); do
		cmd="$project_dir/KernelRidgeRegression/HCP/scripts/CBIG_LiGSR_KRR_workflowHCP.sh -subject_list $subject_list "
		cmd="$cmd -RSFC_file $RSFC_file -y_list $cog_list -covariate_list $covariate_list -FD_file $FD_file -DVARS_file "
		cmd="$cmd $DVARS_file -outdir $outdir -outstem $outstem -seed $seed"
		
		echo $cmd | qsub -V -q circ-spool -l walltime=06:00:00,mem=6GB -m ae -N CBIG_LiGSR_KRR_replication_all_HCP
		
		if [ ! -f $outdir/covariates_${outstem}.mat ] || [ ! -f $outdir/y_${outstem}.mat ]; then
			# wait for the files shared across random splits to be saved
			sleep 3m   
		else
			sleep 3s
		fi
	done
	
	
	##########################
	# 22 personality and task fMRI measures
	##########################
	person_list="$replication_dir/scripts/HCP_lists/Personality_Task_unrestricted.txt"
	covariate_list="$replication_dir/scripts/HCP_lists/covariates_58behaviors.txt"
	outstem=22personality
	
	for seed in $(seq 1 1 20); do
		cmd="$project_dir/KernelRidgeRegression/HCP/scripts/CBIG_LiGSR_KRR_workflowHCP.sh -subject_list $subject_list "
		cmd="$cmd -RSFC_file $RSFC_file -y_list $person_list -covariate_list $covariate_list -FD_file $FD_file -DVARS_file "
		cmd="$cmd $DVARS_file -outdir $outdir -outstem $outstem -seed $seed"
		
		echo $cmd | qsub -V -q circ-spool -l walltime=06:00:00,mem=6GB -m ae -N CBIG_LiGSR_KRR_replication_all_HCP
		
		if [ ! -f $outdir/covariates_${outstem}.mat ] || [ ! -f $outdir/y_${outstem}.mat ]; then
			sleep 3m
		else
			sleep 3s
		fi
	done
	
	
	##########################
	# 23 Social emotional measures
	##########################
	emot_list="$replication_dir/scripts/HCP_lists/Social_Emotion_unrestricted.txt"
	covariate_list="$replication_dir/scripts/HCP_lists/covariates_58behaviors.txt"
	outstem=23emotion
	
	for seed in $(seq 1 1 20); do
		cmd="$project_dir/KernelRidgeRegression/HCP/scripts/CBIG_LiGSR_KRR_workflowHCP.sh -subject_list $subject_list "
		cmd="$cmd -RSFC_file $RSFC_file -y_list $emot_list -covariate_list $covariate_list -FD_file $FD_file -DVARS_file "
		cmd="$cmd $DVARS_file -outdir $outdir -outstem $outstem -seed $seed"
		
		echo $cmd | qsub -V -q circ-spool -l walltime=06:00:00,mem=6GB -m ae -N CBIG_LiGSR_KRR_replication_all_HCP
		
		if [ ! -f $outdir/covariates_${outstem}.mat ] || [ ! -f $outdir/y_${outstem}.mat ]; then
			sleep 3m
		else
			sleep 3s
		fi
	done
	
	
	##########################
	# predict age
	##########################
	age_list="$replication_dir/scripts/HCP_lists/Age_header.txt"
	covariate_list="$replication_dir/scripts/HCP_lists/covariates_age.txt"
	outstem=Age
	
	for seed in $(seq 1 1 20); do
		cmd="$project_dir/KernelRidgeRegression/HCP/scripts/CBIG_LiGSR_KRR_workflowHCP.sh -subject_list $subject_list "
		cmd="$cmd -RSFC_file $RSFC_file -y_list $age_list -covariate_list $covariate_list -FD_file $FD_file -DVARS_file "
		cmd="$cmd $DVARS_file -outdir $outdir -outstem $outstem -seed $seed"
		
		echo $cmd | qsub -V -q circ-spool -l walltime=01:00:00,mem=3GB -m ae -N CBIG_LiGSR_KRR_replication_all_HCP
		
		if [ ! -f $outdir/covariates_${outstem}.mat ] || [ ! -f $outdir/y_${outstem}.mat ]; then
			sleep 3m
		else
			sleep 3s
		fi
	done
	
	
	##########################
	# predict sex
	##########################
	sex_list="$replication_dir/scripts/HCP_lists/Sex_header.txt"
	covariate_list="$replication_dir/scripts/HCP_lists/covariates_sex.txt"
	outstem=Sex
	
	for seed in $(seq 1 1 20); do
		cmd="$project_dir/KernelRidgeRegression/HCP/scripts/CBIG_LiGSR_KRR_workflowHCP.sh -subject_list $subject_list "
		cmd="$cmd -RSFC_file $RSFC_file -y_list $sex_list -covariate_list $covariate_list -FD_file $FD_file -DVARS_file "
		cmd="$cmd $DVARS_file -outdir $outdir -outstem $outstem -seed $seed"
		
		echo $cmd | qsub -V -q circ-spool -l walltime=20:00:00,mem=6GB -m ae -N CBIG_LiGSR_KRR_replication_all_HCP
		
		if [ ! -f $outdir/covariates_${outstem}.mat ] || [ ! -f $outdir/y_${outstem}.mat ]; then
			sleep 3m
		else
			sleep 3s
		fi
	done
	
done


