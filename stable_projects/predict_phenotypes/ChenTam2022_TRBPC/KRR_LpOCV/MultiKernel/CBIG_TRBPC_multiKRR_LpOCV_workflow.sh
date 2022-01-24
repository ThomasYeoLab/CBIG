#!/bin/sh

# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Set default value for parameters in optional input arguments
num_leave_out=3
num_inner_folds=10
group_kernel=none
domain=none
metric=predictive_COD
cluster=CBIG_cluster
# Set default value for other parameters. They are not input arguments as they are very unlikely to be changed, but you
# can still change it here if you want to:
ker_param_file=none
threshold=none

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
echo "feature_files = $feature_files" >> $LF
echo "y_list = $y_list" >> $LF
echo "covariate_list = $covariate_list" >> $LF
echo "FD_file = $FD_file" >> $LF
echo "DVARS_file = $DVARS_file" >> $LF
echo "outdir = $outdir" >> $LF
echo "outstem = $outstem" >> $LF
echo "num_leave_out = $num_leave_out" >> $LF
echo "num_inner_folds = $num_inner_folds" >> $LF
echo "group_kernel = $group_kernel" >> $LF
echo "domain = $domain" >> $LF
echo "metric = $metric" >> $LF

############################
# Call matlab function
############################

# prepare input parameters to the multi-kernel regression

matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; CBIG_TRBPC_multiKRR_LpOCV_prepare_parameters( \
   '$csv_file', '$subject_list', '$feature_files', '$y_list', '$covariate_list', '$FD_file', \
   '$DVARS_file', '$outdir', '$outstem', '$num_leave_out', '$num_inner_folds','$ker_param_file',\
   '$threshold','$group_kernel','$domain','$metric' ); exit; " >> $LF 2>&1

paramfile=${outdir}/param.mat
kernel_folders=${outdir}/kernel_folders.mat
num_test_folds=`cat ${outdir}/num_test_folds.txt`

# run multi-kernel regression
if [ "$cluster" == "none" ];then

	for test_fold in `seq 1 ${num_test_folds}`
	do
		${root_dir}/CBIG_TRBPC_multiKRR_LpOCV_testloop.sh -p ${paramfile} -t ${test_fold} -o ${outdir}
	done

else

	currtime=`date +%s | tail -c 5`
	testLPname=mKRR${currtime}
	for test_fold in `seq 1 ${num_test_folds}`
	do
		cmd="${root_dir}/CBIG_TRBPC_multiKRR_LpOCV_testloop.sh -p ${paramfile} -t ${test_fold} -o ${outdir}"
		errfile=${outdir}/job_err_out/testloop_${test_fold}.err
		outfile=${outdir}/job_err_out/testloop_${test_fold}.out
		${CBIG_CODE_DIR}/setup/CBIG_pbsubmit -walltime 50:00:0 -mem 16gb -joberr ${errfile} -jobout ${outfile}\
		 -cmd "${cmd}" -name ${testLPname}
	done

fi

# wait for all jobs to be done
if [ "$cluster" != "none" ];then
	${CBIG_CODE_DIR}/utilities/scripts/CBIG_check_job_status -n ${testLPname}
fi

# gather results of all cross-validation test folds
matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; load $paramfile; CBIG_TRBPC_multiKRR_LpOCV_pick_optima(\
    ${num_test_folds}, param.outdir, param.outstem); exit; " >> $LF 2>&1
}

#############################
# Function usage
#############################
usage() { echo "
NAME:
	CBIG_TRBPC_multiKRR_LpOCV_workflow.sh

DESCRIPTION:
	This function performs leave-p-out cross-validation workflow for multi-kernel ridge regression.
	
REQUIRED ARGUMENTS:
	-csv_file:   	   csv_file   	    : The csv file containing all behavior and covariants
	-subject_list      subject_list     : Full path of the subject list. Each line of this list is one subject ID.
	-feature_files     feature_files    : Full path of the file storing a a cell array. The cell array contains the
	                                      full path of all the feature files that are required to calculate the 
	                                      kernels. For each feature file, a matrix "feature_mat" is assumed to be
	                                      saved in this file. "feature_mat" can be a 2-D matrix with dimension of
	                                      #features x #subjects or a 3-D matrix with dimension of 
	                                      #ROIs1 x #ROIs2 x #subjects.
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
	-group_kernel      group_kernel     : Full path of the file (.mat) storing how the user would like to group the 
	                                      various kernels for multi KRR. By grouping we mean forcing some kernel to 
	                                      always use the same regularization parameter. Default is no grouping. See 
	                                      help in `CBIG_TRBPC_multiKRR_LpOCV_prepare_parameters.m` for more information
	-domain            domain           : Full path of the file (.mat) storing the domain used to specify the search 
	                                      boundary of the hyperparameter. A 1 x 2 vector "domain" is assumed to be
	                                      saved in this file. The first element of the vector stores the lower bound 
	                                      and the second element stores the upper bound. Default is searching from 0
	                                      to 20.
	-metric            metric           : A string stating which accuracy metric to use when optimize hyperparameters.
	                                      Choose from: corr, COD, predictive_COD, MSE, MSE_norm, MAE, MAE_norm. Default
	                                      is predictive_COD. See help in 
	                                      "CBIG_TRBPC_multiKRR_LpOCV_prepare_parameters.m" for more information.
	-cluster           cluster          : if you do not have a cluster, put it as "none", then the cross-validation will
										  run serially (potentially very slow!). If you have a cluster, don't use this
										  argument. The cross-validation will runs as parallel jobs submitted to your
										  cluster
	
Example:
	$CBIG_CODE_DIR/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/KRR_LpOCV/MultiKernel/
	CBIG_TRBPC_multiKRR_LpOCV_workflow.sh -csv_file xxx/behav_and_covariates.csv -subject_list xxx/subjects_all.txt
	-feature_file xxx/paths_allFC.mat -y_list xxx/all_behav_measures.txt -covariate_list xxx/age_sex_motion.txt
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
		-feature_files)
			feature_files=$1
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
		-group_kernel)
			group_kernel=$1
			shift
			;;
		-domain)
			domain=$1
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
if [ "$feature_files" == "" ]; then
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