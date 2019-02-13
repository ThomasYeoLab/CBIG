#!/bin/sh
#
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

data_csv=/share/users/imganalysis/yeolab/data/GSP_release/scripts/subjects/GSP_extended_140630.csv
ystem=""

root_dir="$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/VarianceComponentModel/scripts/"

main() {
############################
# Echo parameters
############################
mkdir -p $outdir/logs
LF="$outdir/logs/${ystem}.log"
if [ -f $LF ]; then rm $LF; fi

echo "RSFC_file = $RSFC_file" >> $LF
echo "data_csv = $data_csv" >> $LF
echo "trait_list = $trait_list" >> $LF
echo "covariate_list = $covariate_list" >> $LF
echo "FD_file = $FD_file" >> $LF
echo "DVARS_file = $DVARS_file" >> $LF
echo "subject_list = $subject_list" >> $LF
echo "outdir = $outdir" >> $LF
echo "ystem = $ystem" >> $LF
echo "d = $d" >> $LF
echo "num_samples = $num_samples"  >> $LF
echo "rmsub_prefix = $rmsub_prefix" >> $LF

############################
# Call matlab function
############################
matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; CBIG_LiGSR_LME_workflowGSP( \
   '$RSFC_file', '$data_csv', '$trait_list', '$covariate_list', '$FD_file', '$DVARS_file', \
   '$subject_list', '$outdir', '$ystem', '$d', '$num_samples', '$rmsub_prefix' ); exit; " >> $LF 2>&1
   
}


#############################
# Function usage
#############################
usage() { echo "
NAME:
	CBIG_LiGSR_LME_workflowGSP.sh
	
DESCRIPTION:
	This function calls the matlab function '"'CBIG_LiGSR_LME_workflowGSP.m'"' to perform variance component model for the 
	Brain Genomics Superstruct Project (GSP) dataset.
	
REQUIRED ARGUMENTS:
	-RSFC_file       RSFC_file      : The resting-state functional connectivity filename (.mat). It is assumed that a 
	                                  #ROIs x #ROIs x #subjects matrix called '"'corr_mat'"' is saved in this file.
	-trait_list      trait_list     : The full path to a text file of trait names. Each line in this file corresponds to one 
	                                  trait. The trait names should exist as headers in '"'data_csv'"'.
	-covariate_list  covariate_list : The full path to a text file containing the name of all covariates. Each line 
	                                  corresponds to one covariate. Except for the covariates FD and DVARS, all the other 
	                                  covariate names should exist as headers in '"'data_csv'"'.
	-FD_file         FD_file        : The full path to a text file of the mean framewise displacement (FD) of all subjects. 
	                                  The number of lines in '"'FD_file'"' should be the same as the number of lines in 
	                                  '"'sub_list'"'. If the user wants to regress out FD, then the '"'covariate_list'"' 
	                                  should include the string '"'FD'"'. If the covariates do not include FD, this input 
	                                  variable is not needed. The user can pass in NONE.
	-DVARS_file      DVARS_file     : The full path to a text file of the mean DVARS of all subjects. The number of lines in 
	                                  '"'DVARS_file'"' should be the same as the number of lines in '"'subject_list'"'. 
	                                  If the user wants to regress out DVARS, then the '"'covariate_list'"' should include 
	                                  the string '"'DVARS'"' (or '"'DV'"'). If the covariates do not include DVARS, 
	                                  this input variable is not needed. The user can pass in NONE.
	-subject_list    subject_list   : Full path of a text file containing all subject IDs. Each line is one subject ID.
	-outdir          outdir         : Full path of the output directory.
	-d               d              : The number of subjects to be removed for each jackknife sample, e.g. 431.
	-num_samples     num_samples    : The total number of jackknife samples, e.g. 1000.
	-rmsub_prefix    rmsub_prefix   : A string, the prefix to the filename of the subject list to be removed for jackknife 
	                                  samples. The list of removed subject IDs of each jackknife sample will be saved under 
	                                  the output file name 
	                                  \${outdir}/jackknife_lists/\${rmsub_prefix}_choose\${d}_set\${i}.txt,
	                                  where i ranges from  1 to num_samples.
	
OPTIONAL ARGUMENTS:
	-data_csv        data_csv       : The full path of the CSV file containing the covariates and traits of all subjects 
	                                  from the GSP dataset. Default is: 
	                                  /share/users/imganalysis/yeolab/data/GSP_release/scripts/subjects/GSP_extended_140630.csv
	
	-ystem           ystem          : The trait values will be read from '"'data_csv'"' and saved in a .mat file: 
	                                  \${outdir}/y_\${ystem}.mat. For example, if '"'trait_list'"' contains 23 behavioral 
	                                  names, you can set ystem = 23behaviors.

EXAMPLE:
	$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/VarianceComponentModel/scripts/CBIG_LiGSR_LME_workflowGSP.sh 
	-RSFC_file xxx/RSFC_862_Fisher_GSR.mat -trait_list xxx/23behaviors.txt -covariate_list xxx/covariates_23behaviors.txt 
	-FD_file xxx/FD_regressor_862.txt -DVARS_file xxx/DV_regressor_862.txt -subject_list xxx/subject_list_862.txt
	-outdir xxx/ref_output/GSR -ystem 23behaviors -d 431 -num_samples 1000 -rmsub_prefix subjects862

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
		-RSFC_file)
			RSFC_file=$1
			shift;;
		
		-data_csv)  # optional
			data_csv=$1
			shift;;
			
		-trait_list)
			trait_list=$1
			shift;;
		
		-covariate_list)
			covariate_list=$1
			shift;;
		
		-FD_file)
			FD_file=$1
			shift;;
		
		-DVARS_file)
			DVARS_file=$1
			shift;;
		
		-subject_list)
			subject_list=$1
			shift;;
		
		-outdir)
			outdir=$1
			shift;;
		
		-ystem)  #optional
			ystem=$1
			shift;;
		
		-d)
			d=$1
			shift;;
		
		-num_samples)
			num_samples=$1
			shift;;
		
		-rmsub_prefix)
			rmsub_prefix=$1
			shift;;
		
		*) 
			echo "Unknown flag $flag"
			usage; 1>&2; exit 1
			;;
	esac
done


##########################################
# Check Parameters
##########################################
if [ "$RSFC_file" == "" ]; then
	echo "ERROR: RSFC file not specified"
	exit 1
fi
if [ "$data_csv" == "" ]; then
	echo "ERROR: data CSV not specified"
	exit 1
fi
if [ "$trait_list" == "" ]; then
	echo "ERROR: trait list not specified"
	exit 1
fi
if [ "$covariate_list" == "" ]; then
	echo "ERROR: covariate list not specified"
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
if [ "$subject_list" == "" ]; then
	echo "ERROR: subject list not specified"
	exit 1
fi
if [ "$outdir" == "" ]; then
	echo "ERROR: output directory not specified"
	exit 1
fi
if [ "$d" == "" ]; then
	echo "ERROR: number of removed subjects per jackknife sample not specified"
	exit 1
fi
if [ "$num_samples" == "" ]; then
	echo "ERROR: number of jackknife samples not specified"
	exit 1
fi
if [ "$rmsub_prefix" == "" ]; then
	echo "ERROR: prefix for removed subject lists not specified"
	exit 1
fi


####### Execute main()
main


