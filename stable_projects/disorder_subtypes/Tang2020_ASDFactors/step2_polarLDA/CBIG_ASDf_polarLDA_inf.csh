#!/bin/csh

##############################################################
# Inferring factor composition of new participants with polarLDA
##############################################################
# AUTHOR ####################################
# Nanbo Sun 
# 2016/07/27
#############################################
#############################################
# In this script, we: 
# 1) run polarLDA to do inference for new participants
# Example: 
# CBIG_ASDf_polarLDA_inf.csh -corpus doc1.dat -model_dir yourModelDir 
# -factor_num 2 -run_num 1 -output_dir yourOutputDir 
# -output_name infFactorComp -infSettings infSettings.txt -code_dir yourCodeDir
#############################################

# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$Id: CBIG_ASDf_polarLDA_inf.csh, v 1.0 2016/07/27 $'

set PrintHelp = 0;
if($#argv == 0) goto usage_exit;
set n = `echo $argv | grep -e -help | wc -l`
if($n != 0) then
	set PrintHelp = 1;
	goto usage_exit;
endif
set n = `echo $argv | grep -e -version | wc -l`
if($n != 0) then
	echo $VERSION
	exit 0;
endif

set corpus = ""
set model_dir = ""
set factor_num = ""
set run_num = ""
set output_dir = ""
set output_name = ""
set infSettings = ""
set code_dir = ""


goto parse_args;
parse_args_return:

goto check_params;
check_params_return:


##########################################
# Run polarLDA inference 
##########################################	
mkdir -p ${output_dir}

set LF = "${output_dir}/k${factor_num}r${run_num}_${output_name}.log"
echo "See ${LF} for outputs and errors."

touch $LF

date >> $LF
echo Corpus: ${corpus} >> $LF
echo k = ${factor_num} >> $LF
echo Inference settings are: >> $LF
cat ${infSettings} >> $LF

set cmd = "${code_dir}/polarLDA inf \
${infSettings} ${model_dir}/k${factor_num}/r${run_num}/final \
${corpus} ${output_dir}/k${factor_num}r${run_num}_${output_name}"
echo $cmd >> $LF
eval $cmd

exit 1
##########################################
# Parse Arguments 
##########################################	

parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
	set flag = $argv[1]; shift;
	
	switch($flag)
		#corpus name
		case "-corpus":
			if ( $#argv == 0 ) goto arg1err;
			set corpus = $argv[1]; shift;
			breaksw	
			
		#model direcdtory
		case "-model_dir":
			if ( $#argv == 0 ) goto arg1err;
			set model_dir = $argv[1]; shift;
			breaksw
			
		#factor number
		case "-factor_num":
			if ( $#argv == 0 ) goto arg1err;
			set factor_num = $argv[1]; shift;
			breaksw
			
		#initialization number
		case "-run_num"
			if ( $#argv == 0 ) goto arg1err;
			set run_num = $argv[1]; shift;
			breaksw

		#output directory
		case "-output_dir":
			if ( $#argv == 0 ) goto arg1err;
			set output_dir = $argv[1]; shift;
			breaksw
			
		#output name
		case "-output_name":
			if ( $#argv == 0 ) goto arg1err;
			set output_name = $argv[1]; shift;
			breaksw

		#inference settings
		case "-infSettings":
			if ( $#argv == 0 ) goto arg1err;
			set infSettings = $argv[1]; shift;
			breaksw

		#code directory
		case "-code_dir":
			if ( $#argv == 0 ) goto arg1err;
			set code_dir = $argv[1]; shift;
			breaksw
			
		default:
			echo ERROR: Flag $flag unrecognized.
			echo $cmdline
			exit 1
			breaksw
	endsw
end
goto parse_args_return;


##########################################
# Check Parameters
##########################################

check_params:
if ( "$corpus" == "" ) then
	echo "ERROR: corpus not specified"
	exit 1;
endif
if ( "$model_dir" == "" ) then
	echo "ERROR: model directory not specified"
	exit 1;
endif
if ( "$factor_num" == "" ) then
	echo "ERROR: factor number not specified"
	exit 1;
endif
if ( "$run_num" == "" ) then
	echo "ERROR: run number not specified"
	exit 1;
endif
if ( "$output_dir" == "" ) then
	echo "ERROR: output directory not specified"
	exit 1;
endif
if ( "$output_name" == "" ) then
	echo "ERROR: output name not specified"
	exit 1;
endif
if ( "$infSettings" == "" ) then
        echo "ERROR: inference settings file not specified"
        exit 1;
endif
if ( "$code_dir" == "" ) then
        echo "ERROR: code directory not specified"
        exit 1;
endif

goto check_params_return;

##########################################
# ERROR message
##########################################
		
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1
arg2err:
  echo "ERROR: flag $flag requires two arguments"
  exit 1
  
#####################################
# Help
#####################################
usage_exit:
	echo ""
	echo "USAGE: CBIG_polarLDA_inf.csh"
	echo ""
	echo "  -corpus  <corpus>         : corpus name"
	echo "  -model_dir <model_dir>    : model directory (output_dir in CBIG_polarLDA_est.csh)"
	echo "  -factor_num <factor_num>  : number of factors"
	echo "  -run_num  <run_num>       : the run that you want to use to do inference"
	echo "  -output_dir  <output_dir> : output directory"
	echo "  -output_name <output_name>: output file name"
	echo "  -infSettings <infSettings>: inference settings text file"
	echo "  -code_dir <code_dir>      : code directory"
	echo ""
	echo "  -help                     : help"
	echo "  -version                  : version"
	echo ""
	
	if(! $PrintHelp) exit 1;
	
	echo $VERSION
	
	cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'
	
exit 1;
#-------- Everything below is printed as part of help --------#
BEGINHELP

This function will take the output of CBIG_ASDf_polarLDA_est.csh results to do inference for a corpus.
Example: 
CBIG_ASDf_polarLDA_inf.csh -corpus doc1.dat -model_dir yourModelDir -factor_num 2 -run_num 1 
-output_dir yourOutputDir -output_name infFactorComp -infSettings infSettings.txt -code_dir yourCodeDir
