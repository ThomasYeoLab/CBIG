#!/bin/sh
#
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Set default value for parameters in optional input arguments
num_leave_out=3
num_inner_folds=10
lambda_set_file=none
metric=predictive_COD
culster=CBIG_cluster
# Set default value for other parameters. They are not input arguments as they are very unlikely to be changed, but you
# can still change it here if you want to:
ker_param_file=none
threshold_set_file=none

main() {
############################
# Echo parameters
############################
root_dir=`dirname "$(readlink -f "$0")"`
export MATLABPATH=$CBIG_CODE_DIR/setup
project_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
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

############################
# Call matlab function
############################

# prepare input parameters to the single-kernel regression
matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; CBIG_TRBPC_KRR_LpOCV_prepare_parameters( '$csv_file',\
    '$subject_list', '$feature_file', '$y_list', '$covariate_list', '$FD_file', '$DVARS_file', '$outdir', '$outstem', \
	'$num_leave_out','$num_inner_folds', '$ker_param_file','$lambda_set_file','$threshold_set_file',\
	'$metric' ); exit; " >> $LF 2>&1

# run testloop cross validation
paramfile=${outdir}/param.mat
num_test_folds=`cat ${outdir}/num_test_folds.txt`
if [ "$cluster" == "none" ];then
	
	${root_dir}/CBIG_TRBPC_KRR_LpOCV_testloop.sh -p ${paramfile} -n ${num_test_folds} -o ${outdir}

else

	cmd="${root_dir}/CBIG_TRBPC_KRR_LpOCV_testloop.sh -p ${paramfile} -n ${num_test_folds} -o ${outdir}"
	errfile=${outdir}/job_err_out/testloop.err
	outfile=${outdir}/job_err_out/testloop.out
	currtime=`date +%s | tail -c 5`
	testLPname=sKRR${currtime}
	${CBIG_CODE_DIR}/setup/CBIG_pbsubmit -walltime 2:00:0 -mem 8gb -joberr ${errfile} -jobout ${outfile}\
	 -cmd "${cmd}" -name ${testLPname}

fi

# run innerloop cross-validation
if [ "$cluster" == "none" ];then

	for test_fold in `seq 1 ${num_test_folds}`
	do
		${root_dir}/CBIG_TRBPC_KRR_LpOCV_innerloop.sh -p ${paramfile} -t ${test_fold} -o ${outdir}
	done

else
	currtime=`date +%s | tail -c 5`
	innerLPname=sKRR${currtime}
	for test_fold in `seq 1 ${num_test_folds}`
	do
		cmd="${root_dir}/CBIG_TRBPC_KRR_LpOCV_innerloop.sh -p ${paramfile} -t ${test_fold} -o ${outdir}"
		errfile=${outdir}/job_err_out/innerloop_${test_fold}.err
		outfile=${outdir}/job_err_out/innerloop_${test_fold}.out
		${CBIG_CODE_DIR}/setup/CBIG_pbsubmit -walltime 2:00:0 -mem 8gb -joberr ${errfile} -jobout ${outfile}\
		 -cmd "${cmd}" -name ${innerLPname}
	done

fi

# wait for innerloop and testloop to finish
if [ "$cluster" != "none" ];then
	${CBIG_CODE_DIR}/utilities/scripts/CBIG_check_job_status -n ${testLPname}
	${CBIG_CODE_DIR}/utilities/scripts/CBIG_check_job_status -n ${innerLPname}
fi

# gather results of all cross-validation test folds
matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; load $paramfile; CBIG_KRR_pick_optima( \
    param.sub_fold, param.outdir, param.outstem, param.bin_flag, param.ker_param, param.lambda_set,\
    param.threshold_set ); exit; " >> $LF 2>&1
   
}

#############################
# Function usage
#############################
usage() { echo "
NAME:
	CBIG_TRBPC_KRR_LpOCV_workflow.sh

DESCRIPTION:
	This function performs leave-p-out cross-validation workflow for single-kernel ridge regression.
	
REQUIRED ARGUMENTS:
	-csv_file:   	   csv_file   	    : The csv file containing all behavior and covariates
	-subject_list      subject_list     : Full path of the subject list. Each line of this list is one subject ID.
	-feature_file      feature_file     : Full path of the feature file.
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
	-lambda_set_file   lambda_set_file  : Full path of the regularization parameter file (.mat). A vector "lambda_set"
	                                      is assumed to be saved in this file. "lambda_set" is a vector of numbers for
	                                      grid search of lambda (the regularization parameter). If this file is not 
	                                      passed in or set as "none", we will use the default lambda set. Default
										  lambda set is [0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 
										   0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20]
	-metric            metric           : A string stating which accuracy metric to use when optimize hyperparameters.
	                                      Choose from: corr, COD, predictive_COD, MSE, MSE_norm, MAE, MAE_norm. Default
	                                      is predictive_COD. See help in 
	                                      `CBIG_TRBPC_KRR_LpOCV_prepare_parameters.m` for more information
	-cluster           cluster          : if you do not have a cluster, put it as "none", then the cross-validation will
										  run serially (potentially very slow!). If you have a cluster, don't use this
										  argument. The cross-validation will runs as parallel jobs submitted to your
										  cluster
	
Example:
	$CBIG_CODE_DIR/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/KRR_LpOCV/SingleKernel/
	CBIG_TRBPC_KRR_LpOCV_workflow.sh -csv_file xxx/behav_and_covariates.csv -subject_list xxx/subjects_all.txt
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
		-lambda_set_file)
			lambda_set_file=$1
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
	echo "ERROR: feature file not specified"
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
if [ "$lambda_set_file" == "" ]; then
	echo "ERROR: lambda set file not specified"
	exit 1
fi
if [ "$metric" == "" ]; then
	echo "ERROR: optimizing metric not specified"
	exit 1
fi

main



