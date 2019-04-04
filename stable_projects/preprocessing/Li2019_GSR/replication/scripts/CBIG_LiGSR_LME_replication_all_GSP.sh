#!/bin/sh
# This function replicates the variance component model results in the GSP dataset shown in Li et al., 2019
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

test_dir=/mnt/eql/yeo1/CBIG_private_data/replication/stable_projects/preprocessing/Li2019_GSR/VarianceComponentModel/GSP
subject_list="$test_dir/lists/subject_list_862.txt"
FD_file="$test_dir/lists/FD_regressor_862.txt"
DVARS_file="$test_dir/lists/DV_regressor_862.txt"
d=431
num_samples=1000
rmsub_prefix="subjects862"

top_outdir=$1

for pipeline in GSR Baseline ; do
	RSFC_file=$test_dir/RSFC_862_Fisher_${pipeline}.mat
	outdir=$top_outdir/$pipeline
	
	##########################
	# 23 behavioral measures
	##########################
	cog_list="$replication_dir/scripts/GSP_lists/23behaviors.txt"
	covariate_list="$replication_dir/scripts/GSP_lists/covariates_23behaviors.txt"
	ystem=23behaviors
	
	cmd="$project_dir/VarianceComponentModel/scripts/CBIG_LiGSR_LME_workflowGSP.sh -RSFC_file $RSFC_file -trait_list "
	cmd="$cmd $cog_list -covariate_list $covariate_list -FD_file $FD_file -DVARS_file $DVARS_file -subject_list $subject_list"
	cmd="$cmd -outdir $outdir -ystem $ystem -d $d -num_samples $num_samples -rmsub_prefix $rmsub_prefix"
	
	#echo $cmd
	#exit 0
	echo $cmd | qsub -V -q circ-spool -l walltime=01:00:00,mem=4GB,nodes=1:ppn=10 -m ae -N CBIG_LiGSR_LME_replication_all_GSP
	sleep 3s
	
	##########################
	# age
	##########################
	age_list="$replication_dir/scripts/GSP_lists/Age_header.txt"
	covariate_list="$replication_dir/scripts/GSP_lists/covariates_age.txt"
	ystem=Age
	
	cmd="$project_dir/VarianceComponentModel/scripts/CBIG_LiGSR_LME_workflowGSP.sh -RSFC_file $RSFC_file -trait_list "
	cmd="$cmd $age_list -covariate_list $covariate_list -FD_file $FD_file -DVARS_file $DVARS_file -subject_list $subject_list"
	cmd="$cmd -outdir $outdir -ystem $ystem -d $d -num_samples $num_samples -rmsub_prefix $rmsub_prefix"
	
	echo $cmd | qsub -V -q circ-spool -l walltime=00:30:00,mem=4GB,nodes=1:ppn=5 -m ae -N CBIG_LiGSR_LME_replication_all_GSP
	sleep 3s
	
	##########################
	# predict sex
	##########################
	sex_list="$replication_dir/scripts/GSP_lists/Sex_header.txt"
	covariate_list="$replication_dir/scripts/GSP_lists/covariates_sex.txt"
	ystem=Sex
	
	cmd="$project_dir/VarianceComponentModel/scripts/CBIG_LiGSR_LME_workflowGSP.sh -RSFC_file $RSFC_file -trait_list "
	cmd="$cmd $sex_list -covariate_list $covariate_list -FD_file $FD_file -DVARS_file $DVARS_file -subject_list $subject_list"
	cmd="$cmd -outdir $outdir -ystem $ystem -d $d -num_samples $num_samples -rmsub_prefix $rmsub_prefix"
	
	echo $cmd | qsub -V -q circ-spool -l walltime=00:30:00,mem=4GB,nodes=1:ppn=5 -m ae -N CBIG_LiGSR_LME_replication_all_GSP
	sleep 3s
done

