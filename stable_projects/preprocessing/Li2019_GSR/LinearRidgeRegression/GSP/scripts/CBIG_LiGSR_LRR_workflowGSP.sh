#!/bin/sh
#
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

data_csv=$CBIG_LiGSR_REP_GSP_DIR/scripts/subjects/GSP_extended_140630.csv

num_test_folds=20
num_inner_folds=20
evaluations=15
tree=3

root_dir="$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/LinearRidgeRegression/GSP/scripts/"

main() {
############################
# Echo parameters
############################
mkdir -p $outdir/logs/randseed_${seed}
LF="$outdir/logs/randseed_${seed}/${y_name}.log"
if [ -f $LF ]; then rm $LF; fi

echo "data_csv = $data_csv" >> $LF
echo "subject_list = $subject_list" >> $LF
echo "RSFC_file = $RSFC_file" >> $LF
echo "y_name = $y_name" >> $LF
echo "covariate_list = $covariate_list" >> $LF
echo "FD_file = $FD_file" >> $LF
echo "DVARS_file = $DVARS_file" >> $LF
echo "outdir = $outdir" >> $LF
echo "gpso_dir = $gpso_dir" >> $LF
echo "num_test_folds = $num_test_folds" >> $LF
echo "num_inner_folds = $num_inner_folds" >> $LF
echo "seed = $seed" >> $LF
echo "evaluations = $evaluations" >> $LF
echo "tree = $tree" >> $LF

############################
# Call matlab function
############################
matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; CBIG_LiGSR_LRR_workflowGSP( \
   '$data_csv', '$subject_list', '$RSFC_file', '$y_name', '$covariate_list', '$seed', '$num_test_folds', \
   '$num_inner_folds', '$FD_file', '$DVARS_file', '$gpso_dir', '$outdir', '$evaluations', '$tree' ); \
   exit; " >> $LF 2>&1
   
}


#############################
# Function usage
#############################
usage() { echo "
NAME:
    CBIG_LiGSR_KRR_workflowGSP.sh

DESCRIPTION:
    This function calls the matlab function '"'CBIG_LiGSR_KRR_workflowGSP.m'"' to perform linear ridge regression for 
    the Brain Genomics Superstruct Project (GSP) dataset.

REQUIRED ARGUMENTS:
    -subject_list      subject_list     : Full path of the subject list. Each line of this list is one subject ID.
    -RSFC_file         RSFC_file        : Full path of the resting-state functional connectivity matrix.
    -y_name            y_name           : The name of the target measure to be predicted (e.g. behaviors). It should 
                                          correspond to one of the headers stated in \"data_csv\".
    -covariate_list    covariate_list   : Full path of a text list of covariate names (e.g. age, sex, FD) that need to 
                                          be regressed from y (i.e. the measures to be predicted). Each line corresponds
                                          to one covariate name. Except for FD and DVARS, the covariate name should 
                                          correspond to one of the headers stated in \"data_csv\".
    -FD_file           FD_file          : Full path of a text file with the mean FD of each subject. Each line 
                                          corresponds to one subject. The number of lines in FD_file should be the same 
                                          as the number of lines in subject_list.
    -DVARS_file        DVARS_file       : Full path of a text file with the mean DVARS of each subject. Each line 
                                          corresponds to one subject. The number of lines in DVARS_file should be the 
                                          same as the number of lines in subject_list.
    -gpso_dir          gpso_dir         : Full path of the parent directory where the Gaussian process optimization 
                                          package (GPSO; https://github.com/jhadida/gpso) and its dependency Deck 
                                          (https://github.com/jhadida/deck) reposotories stored. 
    -outdir            outdir           : The output directory.
    -seed              seed             : The random seed used for cross-validation fold split.

OPTIONAL ARGUMENTS:
    -data_csv          data_csv         : The CSV file containing the behavioral and demographic information from the 
                                          GSP dataset. If not passed in, the default is 
                                         \${CBIG_LiGSR_REP_GSP_DIR}/scripts/subjects/GSP_extended_140630.csv
    -num_test_folds    num_test_folds   : The number of training-test cross-validation folds. Default is 20.
    -num_inner_folds   num_inner_folds  : The number of inner-loop cross-validation folds split within each training 
                                          fold (for hyperparameter selection). Default is 20.
    -eval              evaluations      : The maximal evaluation times of the objective function. Default is 15. If 
                                          it is too large, the runtime would be too long; if it is too small, the 
                                          objective function may not be able to reach its optimal value. The users need 
                                          to consider this trade-off when setting this variable. For more information, 
                                          please refer to: https://github.com/jhadida/gpso
    -tree              tree             : The depth of the partition tree used to explore hyperparameters. Default is 3.
                                          For more information, please refer to: https://github.com/jhadida/gpso

Example:
    $CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/LinearRidgeRegression/GSP/scripts/\
    CBIG_LiGSR_LRR_workflowGSP.sh 
    -subject_list xxx/subject_list_862.txt -RSFC_file xxx/cort+subcort_new_S1200_953_Fisher.mat -yname 
    Flank_S_CORRpc -covariate_list xxx/covariates_23behaviors.txt -FD_file xxx/FD_regressor_862.txt
    -DVARS_file xxx/DV_regressor_862.txt -outdir xxx/ref_output -seed 1 -gpso_dir ~/code

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
        -data_csv)  # optional
            data_csv=$1
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
        -y_name)
            y_name=$1
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
        -gpso_dir)
            gpso_dir=$1
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
        -eval)
            evaluations=$1
            shift
            ;;
        -tree)
            tree=$1
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

if [ "$data_csv" == "" ]; then
    echo "ERROR: data CSV not specified"
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
if [ "$y_name" == "" ]; then
    echo "ERROR: the name of the measure to be predicted not specified"
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
if [ "$gpso_dir" == "" ]; then
    echo "ERROR: the directory of GPSO respository not specified"
fi
if [ "$evaluations" == "" ]; then
    echo "ERROR: the number of evaluations not specified"
fi
if [ "$tree" == "" ]; then
    echo "ERROR: the depth of tree searching of hyperparameters not specified"
fi




main



