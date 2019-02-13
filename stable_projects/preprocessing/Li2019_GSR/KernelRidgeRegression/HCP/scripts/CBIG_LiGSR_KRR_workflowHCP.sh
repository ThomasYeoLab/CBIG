#!/bin/sh
#
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


HCP_dir="/share/users/imganalysis/yeolab/data/HCP/S1200/scripts"
restricted_csv="$HCP_dir/restricted_hcp_data/RESTRICTED_jingweili_4_12_2017_1200subjects_fill_empty_zygosityGT_by_zygositySR.csv"
unrestricted_csv="$HCP_dir/subject_measures/unrestricted_jingweili_12_7_2017_21_0_16_NEO_A_corrected.csv"
outstem=""

num_test_folds=20
num_inner_folds=20

root_dir="$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/KernelRidgeRegression/HCP/scripts/"


main() {
############################
# Echo parameters
############################
mkdir -p $outdir/logs
LF="$outdir/logs/randseed_${seed}_${outstem}.log"
if [ -f $LF ]; then rm $LF; fi

echo "restricted_csv = $restricted_csv" >> $LF
echo "unrestricted_csv = $unrestricted_csv" >> $LF
echo "subject_list = $subject_list" >> $LF
echo "RSFC_file = $RSFC_file" >> $LF
echo "y_list = $y_list" >> $LF
echo "covariate_list = $covariate_list" >> $LF
echo "FD_file = $FD_file" >> $LF
echo "DVARS_file = $DVARS_file" >> $LF
echo "outdir = $outdir" >> $LF
echo "outstem = $outstem" >> $LF
echo "num_test_folds = $num_test_folds" >> $LF
echo "num_inner_folds = $num_inner_folds" >> $LF
echo "seed = $seed" >> $LF

############################
# Call matlab function
############################
matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; CBIG_LiGSR_KRR_workflowHCP( \
   '$restricted_csv', '$unrestricted_csv', '$subject_list', '$RSFC_file', '$y_list', '$covariate_list', '$FD_file', \
   '$DVARS_file', '$outdir', '$outstem', '$num_test_folds', '$num_inner_folds', '$seed' ); exit; " >> $LF 2>&1
   
}

#############################
# Function usage
#############################
usage() { echo "
NAME:
	CBIG_LiGSR_KRR_workflowHCP.sh

DESCRIPTION:
	This function calls the matlab function '"'CBIG_LiGSR_KRR_workflowHCP.m'"' to perform kernel ridge regression for the 
	Human Connectome Project (HCP) dataset.
	
REQUIRED ARGUMENTS:
	-subject_list      subject_list     : Full path of the subject list. Each line of this list is one subject ID.
	-RSFC_file         RSFC_file        : Full path of the resting-state functional connectivity matrix.
	-y_list            y_list           : A text list of target measures to be predicted (e.g. behaviors). Each line 
	                                      corresponds to one measure name. The name should correspond to one of the headers 
	                                      in either restricted_csv or unrestricted_csv.
	-covariate_list    covariate_list   : A text list of covariate names (e.g. age, sex, FD) that need to be regressed from
	                                      y (i.e. the measures to be predicted). Each line corresponds to one covariate name.
	                                      Except for FD and DVARS, the covariate name stated in this list should correspond 
	                                      to one of the headers in either restricted_csv or unrestricted_csv.
	-FD_file           FD_file          : Full path of a text file with the mean FD of each subject. Each line corresponds to
	                                      one subject. The number of lines in FD_file should be the same as the number of 
	                                      lines in subject_list.
	-DVARS_file        DVARS_file       : Full path of a text file with the mean DVARS of each subject. Each line corresponds
	                                      to one subject. The number of lines in DVARS_file should be the same as the number 
	                                      of lines in subject_list.
	-outdir            outdir           : The output directory.
	-seed              seed             : The random seed used for cross-validation fold split.

OPTIONAL ARGUMENTS:
	-restricted_csv:   restricted_csv   : The restricted CSV file downloaded from the HCP website. If not passed in, the 
	                                      default is 
	                                      '"'/share/users/imganalysis/yeolab/data/HCP/S1200/scripts/restricted_hcp_data/
	                                      RESTRICTED_jingweili_4_12_2017_1200subjects_fill_empty_zygosityGT_by_zygositySR.csv'"'
	-unrestricted_csv  unrestricted_csv : The unrestricted CSV file downloaded from the HCP website. If not passed in, the 
	                                      default is 
	                                      '"'/share/users/imganalysis/yeolab/data/HCP/S1200/scripts/subject_measures/
	                                      unrestricted_jingweili_12_7_2017_21_0_16_NEO_A_corrected.csv'"'
	-outstem           outstem          : A string appended to the filename to specify the output files. For example, if 
	                                      outstem = 58behaviors, then the output files will be named with a suffix of 
	                                      _58behaviors. If not passed in, the default is without any suffix.
	-num_test_folds    num_test_folds   : The number of training-test cross-validation folds. Default is 20.
	-num_inner_folds   num_inner_folds  : The number of inner-loop cross-validation folds split within each training fold 
	                                      (for hyperparameter selection). Default is 20.
	
Example:
	$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/KernelRidgeRegression/HCP/scripts/CBIG_LiGSR_KRR_workflowHCP.sh 
	-subject_list xxx/subject_list_953.txt -RSFC_file xxx/cort+subcort_new_S1200_953_Fisher.mat -y_list 
	xxx/Cognitive_unrestricted.txt -covariate_list xxx/covariates_58behaviors.txt -FD_file xxx/FD_regressor_953.txt
	-DVARS_file xxx/DV_regressor_953.txt -outdir xxx/ref_output -seed 1

Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
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
		-restricted_csv)   # optional
			restricted_csv=$1
			shift
			;;
		-unrestricted_csv)  # optional
			unrestricted_csv=$1
			shift
			;;
		-subject_list)
			subject_list=$1
			shift
			;;
		-RSFC_file)
			RSFC_file=$1
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
		-num_test_folds)
			num_test_folds=$1
			shift
			;;
		-num_inner_folds)
			num_inner_folds=$1
			shift
			;;
		-seed)
			seed=$1
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

if [ "$restricted_csv" == "" ]; then
	echo "ERROR: restricted CSV not specified"
	exit 1
fi
if [ "$unrestricted_csv" == "" ]; then
	echo "ERROR: unrestricted CSV not specified"
	exit 1
fi
if [ "$subject_list" == "" ]; then
	echo "ERROR: subject list not specified"
	exit 1
fi
if [ "$RSFC_file" == "" ]; then
	echo "ERROR: RSFC file not specified"
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
if [ "$num_test_folds" == "" ]; then
	echo "ERROR: the number of test folds not specified"
	exit 1
fi
if [ "$num_inner_folds" == "" ]; then
	echo "ERROR: the number of inner-loop folds not specified"
	exit 1
fi
if [ "$seed" == "" ]; then
	echo "ERROR: the random seed not specified"
	exit 1
fi




main



