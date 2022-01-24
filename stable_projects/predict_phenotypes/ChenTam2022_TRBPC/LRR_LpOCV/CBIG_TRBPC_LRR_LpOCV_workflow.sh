#!/bin/sh
#
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Set default value for parameters in optional input arguments
num_leave_out=3
num_inner_folds=10
domain=none
eval=none
tree=none
metric=predictive_COD
fold_number=20
cluster=CBIG_cluster

main() {
############################
# Echo parameters
############################
root_dir=`dirname "$(readlink -f "$0")"`
export MATLABPATH=$CBIG_CODE_DIR/setup
mkdir -p $outdir/logs
LF="$outdir/logs/${outstem}.log"
if [ -f $LF ]; then rm $LF; fi

echo "csv_file = $csv_file" >> $LF
echo "subject_list = $subject_list" >> $LF
echo "feature_file = $feature_file" >> $LF
echo "y_list = $y_list" >> $LF
echo "covariate_list = $covariate_list" >> $LF
echo "FD_file = $FD_file" >> $LF
echo "DVARS_file = $DVARS_file" >> $LF
echo "outdir = $outdir" >> $LF
echo "outstem = $outstem" >> $LF
echo "num_leave_out = $num_leave_out" >> $LF
echo "num_inner_folds = $num_inner_folds" >> $LF
echo "domain = $domain" >> $LF
echo "eval = $eval" >> $LF
echo "tree = $tree" >> $LF
echo "metric = $metric" >> $LF

############################
# Call matlab function
############################

# prepare input parameters to the linear ridge regression

matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; CBIG_TRBPC_LRR_LpOCV_prepare_parameters( \
   '$csv_file', '$subject_list', '$feature_file', '$y_list', '$covariate_list', '$FD_file', \
   '$DVARS_file', '$outdir', '$outstem', '$num_leave_out', '$num_inner_folds','$domain','$eval','$tree',\
   '$metric' ); exit; " >> $LF 2>&1

paramfile=${outdir}/param.mat
num_test_folds=`cat ${outdir}/num_test_folds.txt`
num_score=`cat $y_list | wc -l`

# run LRR
if [ "$cluster" == "none" ];then

	for score_ind in `seq 1 ${num_score}`
	do
		for fold_start in `seq 1 ${fold_number} ${num_test_folds}`
		do
			${root_dir}/CBIG_TRBPC_LRR_LpOCV_job.sh -o ${outdir} -b ${score_ind} -s ${fold_start}\
			-n ${fold_number}
		done
	done

else
	currtime=`date +%s | tail -c 5`
	LRRname=LRR${currtime}
	for score_ind in `seq 1 ${num_score}`
	do
		for fold_start in `seq 1 ${fold_number} ${num_test_folds}`
		do
			cmd="${root_dir}/CBIG_TRBPC_LRR_LpOCV_job.sh"
			cmd="$cmd -o ${outdir} -b ${score_ind} -s ${fold_start} -n ${fold_number}"
			errfile=${outdir}/job_err_out/LRR_behav${score_ind}_foldStart${fold_start}.err
			outfile=${outdir}/job_err_out/LRR_behav${score_ind}_foldStart${fold_start}.err
			${CBIG_CODE_DIR}/setup/CBIG_pbsubmit -walltime 6:00:0 -mem 24gb -joberr ${errfile} -jobout ${outfile}\
			 -cmd "${cmd}" -name ${LRRname}
		done
	done
fi
# wait for all jobs to be done
if [ "$cluster" != "none" ];then
	${CBIG_CODE_DIR}/utilities/scripts/CBIG_check_job_status -n ${LRRname}
fi

# gather results of all cross-validation test folds
matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; load $paramfile; CBIG_TRBPC_LRR_LpOCV_pick_optima(\
    ${num_test_folds},${num_score},param.outdir, param.outstem); exit; " >> $LF 2>&1
}

#############################
# Function usage
#############################
usage() { echo "
NAME:
	CBIG_TRBPC_LRR_LpOCV_workflow.sh

DESCRIPTION:
	This function performs leave-p-out cross-validation workflow for linear ridge regression.
	
REQUIRED ARGUMENTS:
	-csv_file:   	   csv_file   	    : The csv file containing all behavior and covariates
	-subject_list      subject_list     : Full path of the subject list. Each line of this list is one subject ID.
	-feature_file      feature_file     : Full path of the feature files
	-y_list            y_list           : A text list of target measures to be predicted (e.g. behaviors). Each line 
	                                      corresponds to one measure name. The name should correspond to one of the 
	                                      headers in the csv_file.
	-covariate_list    covariate_list   : A text list of covariate names (e.g. age, sex, FD) that need to be regressed
	                                      from y (i.e. the measures to be predicted). Each line corresponds to one 
	                                      covariate name. The covariate name stated in this list should 
	                                      correspond to one of the headers in the csv_file.
	-FD_file           FD_file          : Full path of a text file with the mean FD of each subject. Each line 
	                                      corresponds to one subject in the same order as subject_list. Used if FD is 
										  in covariate_list but not a column in csv_file. Put "none" if FD is already 
										  in the csv_file.
	-DVARS_file        DVARS_file       : Full path of a text file with the mean DVARS of each subject. Each line
	                                      corresponds to one subject in the same order as subject_list. Used if DVARS 
										  is in covariate_list but not a column in csv_file. Put "none" if DVARS is 
										  already in the csv_file.
	-outdir            outdir           : The output directory.

OPTIONAL ARGUMENTS:
	-outstem           outstem          : A string appended to the filename to specify the output files. For example, 
	                                      if outstem = all_score, then the output files will be named with a suffix 
	                                      of all_score. If not passed in, the default is without any suffix.
	-num_leave_out     num_leave_out    : The value of p in the leave-p-out cross validation. To choose p note that:
	                                      1): if p is too small then we don't have enough number of cross-validation 
	                                      folds.
	                                      2): if p is too large so we have enough number of subjects in the training 
	                                      set.
	                                      Default is 3 for the ABCD task-rest prediction project.
	-num_inner_folds   num_inner_folds  : The number of inner-loop cross-validation folds split within each training 
	                                      fold (for hyperparameter selection). Default is 10.
	-eval              eval             : The maximal evaluation times of the objective function (a scalar). If it is 
	                                      too large, the runtime would be too long; if it is too small, the objective 
	                                      function may not be able to reach its optimal value. For more information, 
	                                      please refer to: https://github.com/jhadida/gpso Default is 15
	-domain            domain           : The search domain of parameters used by the Gaussian process optimization 
	                                      algorithm (dimension: #hyperparameters x 2). Default is [0.001 0.1; 3 8], 
	                                      where 0.001-0.1 is the search domain for the feature selection threshold,
	                                      and 3-8 is the search domain for the L2 regularization hyperparameter
	                                      (after taking logarithm). For more information,
	                                      please refer to: https://github.com/jhadida/gpso
	-tree              tree             : A scalar, the depth of the partition tree used to explore hyperparameters. 
	                                      Default is 3. For more information, please refer to:
	                                      https://github.com/jhadida/gpso
	-metric            metric           : A string stating which accuracy metric to use when optimize hyperparameters.
	                                      Choose from: corr, COD, predictive_COD, MSE, MSE_norm, MAE, MAE_norm. Default
	                                      is predictive_COD. See help in 
	                                      `CBIG_TRBPC_multiKRR_LpOCV_prepare_parameters.m` for more information
	-cluster           cluster          : if you do not have a cluster, put it as "none", then the cross-validation will
										  run serially (potentially very slow!). If you have a cluster, don't use this
										  argument. The cross-validation will runs as parallel jobs submitted to your
										  cluster
	
Example:
	$CBIG_CODE_DIR/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/LRR_LpOCV/CBIG_TRBPC_LRR_LpOCV_workflow.sh 
	-csv_file xxx/behav_and_covariates.csv -subject_list xxx/subjects_all.txt
	-feature_file xxx/RSFC_all_subjects.mat -y_list xxx/all_behav_measures.txt -covariate_list xxx/age_sex_motion.txt
	-FD_file none -DVARS_file none -outdir xxx/ref_output

" 1>&2; exit 1; }

##########################################
# ERROR message
##########################################	
arg1err() {
	echo "ERROR: flag $flag requires one argument"
	exit 1
}

##########################################
# Parse Arguments 
##########################################
# Display help message if no argument is supplied
if [ $# -eq 0 ]; then
	usage; 1>&2; exit 1
fi

while [[ $# -gt 0 ]]; do
	flag=$1; shift;
	
	case $flag in
		-csv_file)   # optional
			csv_file=$1
			shift
			;;
		-subject_list)
			subject_list=$1
			shift
			;;
		-feature_file)
			feature_file=$1
			shift
			;;
		-y_list)
			y_list=$1
			shift
			;;
		-covariate_list)
			covariate_list=$1
			shift
			;;
		-FD_file)
			FD_file=$1
			shift
			;;
		-DVARS_file)
			DVARS_file=$1
			shift
			;;
		-outdir)
			outdir=$1
			shift
			;;
		-outstem)
			outstem=$1
			shift
			;;
		-num_leave_out)
			num_leave_out=$1
			shift
			;;
		-num_inner_folds)
			num_inner_folds=$1
			shift
			;;
		-domain)
			domain=$1
			shift
			;;
		-eval)
			eval=$1
			shift
			;;
		-tree)
			tree=$1
			shift
			;;
		-metric)
			metric=$1
			shift
			;;
		-cluster)
			cluster=$1
			shift
			;;
		*) 
			echo "Unknown flag $flag"
			usage; 1>&2; exit 1
			;;
	esac
done

##########################################
# Check Parameters
##########################################

if [ "$csv_file" == "" ]; then
	echo "ERROR: restricted CSV not specified"
	exit 1
fi
if [ "$subject_list" == "" ]; then
	echo "ERROR: subject list not specified"
	exit 1
fi
if [ "$feature_file" == "" ]; then
	echo "ERROR: feature files not specified"
	exit 1
fi
if [ "$y_list" == "" ]; then
	echo "ERROR: the list of measures to be predicted not specified"
	exit 1
fi
if [ "$covariate_list" == "" ]; then
	echo "ERROR: the list of covariate names not specified"
	exit 1
fi
if [ "$FD_file" == "" ]; then
	echo "ERROR: FD file not specified"
	exit 1
fi
if [ "$DVARS_file" == "" ]; then
	echo "ERROR: DVARS file not specified"
	exit 1
fi
if [ "$outdir" == "" ]; then
	echo "ERROR: output directory not specified"
	exit 1
fi
if [ "$outstem" == "" ]; then
	echo "ERROR: output stem not specified"
	exit 1
fi
if [ "$num_leave_out" == "" ]; then
	echo "ERROR: the number of leave out sites not specified"
	exit 1
fi
if [ "$num_inner_folds" == "" ]; then
	echo "ERROR: the number of inner-loop folds not specified"
	exit 1
fi

main
