#!/bin/sh
#
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

HCP_dir="/share/users/imganalysis/yeolab/data/HCP/S1200/scripts"
restricted_csv="$HCP_dir/restricted_hcp_data/RESTRICTED_jingweili_4_12_2017_1200subjects_fill_empty_zygosityGT_by_zygositySR.csv"
unrestricted_csv="$HCP_dir/subject_measures/unrestricted_jingweili_12_7_2017_21_0_16_NEO_A_corrected.csv"
unrelated=1
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
echo "restricted_csv = $restricted_csv" >> $LF
echo "unrestricted_csv = $unrestricted_csv" >> $LF
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
matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; CBIG_LiGSR_LME_workflowHCP( \
   '$RSFC_file', '$restricted_csv', '$unrestricted_csv', '$trait_list', '$covariate_list', '$FD_file', '$DVARS_file', \
   '$subject_list', '$unrelated', '$outdir', '$ystem', '$d', '$num_samples', '$rmsub_prefix' ); exit; " >> $LF 2>&1
   
}


#############################
# Function usage
#############################
usage() { echo "
NAME:
	CBIG_LiGSR_LME_workflowHCP.sh
	
DESCRIPTION:
	This function calls the matlab function '"'CBIG_LiGSR_LME_workflowHCP.m'"' to perform variance component model for the 
	Human Connectome Project (HCP) dataset.
	
REQUIRED ARGUMENTS:
	-RSFC_file         RSFC_file        : The resting-state functional connectivity filename (.mat). It is assumed that a 
	                                      #ROIs x #ROIs x #subjects matrix called '"'corr_mat'"' is saved in this file.
	-trait_list        trait_list       : The full path to a text file of trait names. Each line in this file corresponds to 
	                                      one trait. The trait names should exist as headers in '"'data_csv'"'.
	-covariate_list    covariate_list   : The full path to a text file containing the name of all covariates. Each line 
	                                      corresponds to one covariate. Except for the covariates FD and DVARS, all the other 
	                                      covariate names should exist as headers in '"'data_csv'"'.
	-FD_file           FD_file          : The full path of the mean framewise displacement (FD) of all subjects. The number 
	                                      of lines in '"'FD_file'"' should be the same as the number of lines in '"'sub_list'"'. 
	                                      If the user wants to regress out FD, then '"'FD'"' should be included in 
	                                      '"'covariate_list'"'. If the covariates do not include FD, this input variable is 
	                                      not needed. The user can pass in NONE.
	-DVARS_file        DVARS_file       : The full path of the mean DVARS of all subjects. The number of lines in '"'DVARS_file'"' 
	                                      should be the same as the number of lines in '"'subject_list'"'. If the user wants 
	                                      to regress out '"'DVARS'"', then '"'DVARS'"' (or '"'DV'"') should be included in 
	                                      '"'covariate_list'"'. If the covariates do not include DVARS, this input variable 
	                                      is not needed. The user can pass in NONE.
	-subject_list      subject_list     : Full path of a text file containing all subjects. Each line is one subject ID.
	-outdir            outdir           : Full path of the output directory.
	-d                 d                : The number of subjects to be removed for each jackknife sample, e.g. 209.
	-num_samples       num_samples      : The total number of jackknife samples, e.g. 1000.
	-rmsub_prefix      rmsub_prefix     : A string. The prefix to the filename of the subject list to be removed for 
	                                      jackknife samples. The list of removed subject IDs of each jackknife sample will 
	                                      be saved under the output file name 
	                                      \${outdir}/jackknife_lists/\${rmsub_prefix}_choose\${d}_set\${i}.txt,
	                                      where i ranges from  1 to num_samples.
	
OPTIONAL ARGUMENTS:
	-restricted_csv    restricted_csv   : The restricted CSV file downloaded from the HCP website. If not passed in, the 
	                                      default is 
	                                      '"'/share/users/imganalysis/yeolab/data/HCP/S1200/scripts/restricted_hcp_data/
	                                      RESTRICTED_jingweili_4_12_2017_1200subjects_fill_empty_zygosityGT_by_zygositySR.csv'"'
	-unrestricted_csv  unrestricted_csv : The unrestricted CSV file downloaded from the HCP website. If not passed in, the 
	                                      default is 
	                                      '"'/share/users/imganalysis/yeolab/data/HCP/S1200/scripts/subject_measures/
	                                      unrestricted_jingweili_12_7_2017_21_0_16_NEO_A_corrected.csv'"'
	-unrelated         unrelated        : 0 or 1. 1 means all subjects in the \$subject_list are unrelated. 0 means there is 
	                                      a family structure within the subjects in \$subject_list, and the family 
	                                      information will be read from \$restricted_csv. Default is 1.
	-ystem             ystem            : The trait values will be read from '"'restricted_csv'"' and '"'unrestricted_csv'"',
	                                      and saved in a .mat file: 
	                                      \${outdir}/y_\${ystem}.mat. For example, if '"'trait_list'"' contains 13 
	                                      cognitive behavioral names, you can set ystem = 13cognitive.

EXAMPLE:
	$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/VarianceComponentModel/scripts/CBIG_LiGSR_LME_workflowHCP.sh 
	-RSFC_file xxx/RSFC_953_unrelated_419_Fisher_GSR.mat -trait_list xxx/13cognitive.txt -covariate_list 
	xxx/covariates_58behaviors.txt -FD_file xxx/FD_regressor_953_unrelated_419.txt -DVARS_file 
	xxx/DV_regressor_953_unrelated_419.txt -subject_list xxx/subject_list_953_unrelated_419.txt -outdir xxx/ref_output/GSR 
	-ystem 13cognitive -d 208 -num_samples 1000 -rmsub_prefix subjects953_unrelated419

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
		
		-restricted_csv)  # optional
			restricted_csv=$1
			shift;;
			
		-unrestricted_csv)  # optional
			unrestricted_csv=$1
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
			
		-unrelated)
			unrelated=$1
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
if [ "$restricted_csv" == "" ]; then
	echo "ERROR: restricted CSV not specified"
	exit 1
fi
if [ "$unrestricted_csv" == "" ]; then
	echo "ERROR: unrestricted CSV not specified"
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
if [ "$unrelated" == "" ]; then
	echo "ERROR: whether the input subjects are related or unrelated is not specified"
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


