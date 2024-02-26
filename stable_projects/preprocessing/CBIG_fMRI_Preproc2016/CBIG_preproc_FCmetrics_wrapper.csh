#!/bin/csh -f

###################################################################
# Compute FC metrics for each subject
###################################################################
# Author ##########################################################
# Jingwei Li
# Oct. 14, 2017
###################################################################
# This function will:
# 1. Extract the 19 subcortical ROIs from aseg in subject-specific 
#    functional space.
# 2. Generate input lists for the Matlab function
# 3. Call the Matlab function to compute FC metrics of cortical 
#    and 19 subcortical ROIs. Currently we only support to compute
#    static Pearson's correlation by using "-Pearson_r" option.
#    In the future we will include other types of metrics.
###################################################################
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$Id: CBIG_preproc_FCmetrics_wrapper.csh, v 1.0 2017/10/14 $'

set n = `echo $argv | grep -e -help | wc -l`

# if there is -help option 
if( $n != 0 ) then
    echo $VERSION
    # print help
    cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'
    exit 0;
endif

# if there is no arguments
if( $#argv == 0 ) then
    echo $VERSION
    # print help
    cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'
    echo "WARNING: No input arguments. See above for a list of available input arguments."
    exit 0;
endif


set n = `echo $argv | grep -e -version | wc -l`
if($n != 0) then
    echo $VERSION
    exit 0;
endif


set subject = ""
set sub_dir = ""
set bold = ""
set BOLD_stem = ""
set surf_stem = ""
set outlier_stem = ""
set lh_cortical_ROIs_file = ""
set rh_cortical_ROIs_file = ""
set parcellation_name = ""
set censor = 0
set Pearson_r = 0
set nocleanup = 0; # Default clean up intermediate file
set parcellation_type = ""
set network = ""
set res = ""
set matlab_runtime = 0 # Default not running on MATLAB Runtime

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

set root_dir = `python -c "import os; print(os.path.realpath('$0'))"`
set root_dir = `dirname $root_dir`


###############################
# check if matlab exists
###############################
if ( $matlab_runtime == 0 ) then
    set MATLAB=`which $CBIG_MATLAB_DIR/bin/matlab`
    if ($status) then
        echo "ERROR: could not find MATLAB"
        exit 1;
    endif
else
    set MATLAB=$MATLAB_RUNTIME_DIR
    if ($status) then
        echo "ERROR: could not find MATLAB Runtime"
        exit 1;
    endif
    echo "Setting up environment variables for MATLAB Runtime"
    setenv LD_LIBRARY_PATH ${MATLAB}/runtime/glnxa64:${MATLAB}/bin/glnxa64:${MATLAB}/sys/os/glnxa64:${MATLAB}/sys/opengl/lib/glnxa64:${LD_LIBRARY_PATH}
    # check if MATLAB Runtime utilities folder exists
    set matlab_runtime_util = "${root_dir}/matlab_runtime/utilities"
    if ( ! -d $matlab_runtime_util ) then
        echo "ERROR: MATLAB Runtime utilities folder does not exist!"
        exit 1;
    endif
endif

cd $sub_dir/$subject


###############################
# create log file
###############################
if (! -e logs) then
     mkdir -p logs
endif
set LF = $sub_dir/$subject/logs/CBIG_preproc_FCmetrics_wrapper.log
set all_LF = $sub_dir/$subject/logs/CBIG_preproc_fMRI_preprocess.log
set FC_done = `grep -q "\[CBIG_preproc_FC_metrics\]: Done!" $all_LF;echo $?`

if ( "$FC_done" == 0 ) then
    echo "Computing another FC metrics with different config" >> $LF
else
    if( -e $LF ) then
        rm $LF
    endif
    touch $LF
endif

echo "[FC metrics]: logfile = $LF"
echo "Compute FC metrics" >> $LF
echo "[CMD]: CBIG_preproc_FCmetrics_wrapper.csh $cmdline"   >>$LF


###############################
# set parcellation name and ROIs file
###############################
if ( "${parcellation_name}" == "" ) then
    if (( "$parcellation_type" == "" ) || ( "$network" == "" ) || ( "$res" == "" )) then 
        echo "[WARNING]: Default parcellation has been changed to Yan2023 400Parcels with Kong17 networks" \
"from Schaefer2018 400Parcels with Yeo17Networks. Please check!"  |& tee -a $LF
    endif
    if ( "$parcellation_type" == "" ) then
        set parcellation_type = 'Yan'
        echo "Using Yan2023 as default parcellation"  |& tee -a $LF
    endif
    if ( "$network" == "" ) then
        set network = 'Kong17'
        echo "Using Kong2022 17networks as default network assignment"  |& tee -a $LF
    endif
    if ( "$res" == "" ) then
        set res = '400'
        echo "Using 400 Parcels as default resolution"  |& tee -a $LF
    endif

    if ( "$parcellation_type" == 'Yan' ) then
        if ( "${network}" == "Yeo7" ) then
            set network_name = 'Yeo2011_7Networks'
        else if ( "${network}" == "Yeo17" ) then
            set network_name = 'Yeo2011_17Networks'
        else if ( "${network}" == "Kong17" ) then
            set network_name = 'Kong2022_17Networks'
        endif
        set network_folder = `echo "${network}" | tr "[:upper:]" "[:lower:]"`
        set lh_cortical_ROIs_file = "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Yan2023_homotopic/parcellations/FreeSurfer/fsaverage6/label"
        set rh_cortical_ROIs_file = "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Yan2023_homotopic/parcellations/FreeSurfer/fsaverage6/label"
        set lh_cortical_ROIs_file = "${lh_cortical_ROIs_file}/${network_folder}/lh.${res}Parcels_${network_name}.annot"
        set rh_cortical_ROIs_file = "${rh_cortical_ROIs_file}/${network_folder}/rh.${res}Parcels_${network_name}.annot"
    else if ( "$parcellation_type" == 'Schaefer' ) then
        if ( "${network}" == "Yeo7" ) then
            set network_name = '7Networks_order'
        else if ( "${network}" == "Yeo17" ) then 
            set network_name = '17Networks_order'
        else if ( "${network}" == "Kong17" ) then 
            set network_name = 'Kong2022_17Networks_order'
        endif
        set lh_cortical_ROIs_file = "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/"
        set lh_cortical_ROIs_file = "${lh_cortical_ROIs_file}FreeSurfer5.3/fsaverage6/label"
        set rh_cortical_ROIs_file = "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/"
        set rh_cortical_ROIs_file = "${rh_cortical_ROIs_file}FreeSurfer5.3/fsaverage6/label"
        set lh_cortical_ROIs_file = "${lh_cortical_ROIs_file}/lh.Schaefer2018_${res}Parcels_${network_name}.annot"
        set rh_cortical_ROIs_file = "${rh_cortical_ROIs_file}/rh.Schaefer2018_${res}Parcels_${network_name}.annot"
    endif
    set parcellation_name = "${parcellation_type}_${res}Parcels_${network}Networks_order"
endif


###############################
# set folders
###############################
set boldfolder = "$sub_dir/$subject/bold"
echo "[FC metrics]: boldfolder = $boldfolder" |& tee -a $LF
set maskfolder = "$boldfolder/mask"
echo "[FC metrics]: maskfolder = $maskfolder" |& tee -a $LF
set surffolder = "$sub_dir/$subject/surf"
echo "[FC metrics]: surffolder = $surffolder" |& tee -a $LF
set qcfolder = "$sub_dir/$subject/qc"
echo "[FC metrics]: qcfolder = $qcfolder" |& tee -a $LF
set FCmetrics_folder = "$sub_dir/$subject/FC_metrics"
echo "[FC metrics]: FCmetrics_folder = $FCmetrics_folder" |& tee -a $LF


#########################################
# Extract subcortical labels from aseg in functional space
#########################################
echo "======================= Extract 19 subcortical labels from aseg in functional space ========================" \
|& tee -a $LF

set aseg_func = "$maskfolder/$subject.func.aseg.nii"
mkdir -p $FCmetrics_folder/ROIs
set subcortex_func_vol_bin = "$FCmetrics_folder/ROIs/${subject}.subcortex.19aseg.func.bin.nii.gz"
set subcortex_func_vol = "$FCmetrics_folder/ROIs/${subject}.subcortex.19aseg.func.nii.gz"

if ( ! -e $subcortex_func_vol_bin ) then
    set cmd = "mri_binarize --i $aseg_func --match 8 10 11 12 13 16 17 18 26 28 47 49 50 51 52 53 54 58 60 --o \
$subcortex_func_vol_bin"
    echo $cmd |& tee -a $LF
    eval $cmd |& tee -a $LF

    if ( ! -e $subcortex_func_vol_bin ) then
        echo "ERROR: Binarizing subcortical segments in functional aseg failed." |& tee -a $LF
        exit 1
    else
        echo "[FC metrics]: $subcortex_func_vol_bin is created." |& tee -a $LF
    endif
else
    echo "[FC metrics]: $subcortex_func_vol_bin already exists." |& tee -a $LF
endif

if ( ! -e $subcortex_func_vol ) then
    set cmd = "fslmaths $aseg_func -mul $subcortex_func_vol_bin $subcortex_func_vol"
    echo $cmd |& tee -a $LF
    eval $cmd |& tee -a $LF

    if ( ! -e $subcortex_func_vol_bin ) then
        echo "ERROR: Converting binarized subcortical volume to subcortical label volume failed." |& tee -a $LF
        exit 1
    else
        echo "[FC metrics]: $subcortex_func_vol is created." |& tee -a $LF
    endif
else
    echo "[FC metrics]: $subcortex_func_vol already exists." |& tee -a $LF
endif

echo "=================== Extracting 19 subcortical labels from aseg in functional space finished ==================" \
|& tee -a $LF


################################################
# Generate input lists for matlab function
################################################
echo "============================= Generate inout lists for matlab function ===========================" |& tee -a $LF
mkdir -p $FCmetrics_folder/lists
set lh_surf_data_list = "$FCmetrics_folder/lists/lh.${subject}${surf_stem}_list.txt"
set rh_surf_data_list = "$FCmetrics_folder/lists/rh.${subject}${surf_stem}_list.txt"
set vol_data_list = "$FCmetrics_folder/lists/${subject}${BOLD_stem}_list.txt"
if ( $censor == 1 ) then
    set discard_frames_list = "$FCmetrics_folder/lists/${subject}${outlier_stem}_list.txt"
else
    set discard_frames_list = "NONE"
endif

if ( -e $lh_surf_data_list ) then
    rm $lh_surf_data_list
endif
if ( -e $rh_surf_data_list ) then
    rm $rh_surf_data_list
endif
if ( -e $vol_data_list ) then
    rm $vol_data_list
endif
if ( $censor == 1 ) then
    if ( -e $discard_frames_list ) then
        rm $discard_frames_list
    endif
endif

set lh_tmp = ""
set rh_tmp = ""
set vol_tmp = ""
set discard_tmp = ""
foreach run ($bold)
    set lh_tmp = "$lh_tmp $surffolder/lh.${subject}_bld${run}${surf_stem}.nii.gz"
    set rh_tmp = "$rh_tmp $surffolder/rh.${subject}_bld${run}${surf_stem}.nii.gz"
    set vol_tmp = "$vol_tmp $boldfolder/$run/${subject}_bld${run}${BOLD_stem}.nii.gz"
    if ( $censor == 1 ) then
        set discard_tmp = "$discard_tmp $qcfolder/${subject}_bld${run}${outlier_stem}"
    endif
end
echo $lh_tmp > $lh_surf_data_list
echo $rh_tmp > $rh_surf_data_list
echo $vol_tmp > $vol_data_list
if ( $censor == 1 ) then
    echo $discard_tmp > $discard_frames_list
endif

if ( $censor == 1 ) then
    if ( -e $lh_surf_data_list && -e $rh_surf_data_list && -e $vol_data_list && -e $discard_frames_list ) then
        echo "[FC metrics]: Generating input lists for matlab function sucessed." |& tee -a $LF
    else
        echo "ERROR: Generating input lists for matlab function failed." |& tee -a $LF
        exit 1
    endif
else
    if ( -e $lh_surf_data_list && -e $rh_surf_data_list && -e $vol_data_list ) then
        echo "[FC metrics]: Generating input lists for matlab function sucessed." |& tee -a $LF
    else
        echo "ERROR: Generating input lists for matlab function failed." |& tee -a $LF
        exit 1
    endif
endif
echo "========================= Generate input lists for matlab function finished ======================" |& tee -a $LF


#######################################
# Compute FC metrics
#######################################
echo "======================================= Compute FC metrics =======================================" |& tee -a $LF
if ( $Pearson_r == 1 ) then
    set output_dir = "$FCmetrics_folder/Pearson_r"
    set output_prefix = "${subject}${surf_stem}"

    # check if final output exist
    if ( -e $output_dir/${output_prefix}_all2all_${parcellation_name}_with_19Subcortical.mat ) then
        echo "[FC metrics]: The output file $output_dir/${output_prefix}_all2all.mat already exists. Skip ..." |& tee -a $LF
    else
        set matlab_args = "'$lh_cortical_ROIs_file' '$rh_cortical_ROIs_file' '$subcortex_func_vol'"
        set matlab_args = "${matlab_args} '$lh_surf_data_list' '$rh_surf_data_list' '$vol_data_list'"
        set matlab_args = "${matlab_args} '$discard_frames_list' 'Pearson_r' '$output_dir' '$output_prefix'"
        if ( $matlab_runtime == 0 ) then
            set cmd = ( $MATLAB -nodesktop -nodisplay -nosplash -r '"' 'addpath(genpath('"'"${root_dir}'/utilities'"'"'))'; )
            set cmd = ( $cmd CBIG_preproc_FCmetrics $matlab_args; exit; '"' );
        else
            set cmd = ( ${matlab_runtime_util}/CBIG_preproc_FCmetrics $matlab_args )
        endif
        echo $cmd |& tee -a $LF
        eval $cmd |& tee -a $LF

        if ( $nocleanup == 0 ) then
            rm $output_dir/${output_prefix}_lh2lh.mat
            rm $output_dir/${output_prefix}_lh2rh.mat
            rm $output_dir/${output_prefix}_rh2rh.mat
            rm $output_dir/${output_prefix}_lh2subcort.mat
            rm $output_dir/${output_prefix}_rh2subcort.mat
            rm $output_dir/${output_prefix}_subcort2subcort.mat
            mv $output_dir/${output_prefix}_all2all.mat $output_dir/${output_prefix}_all2all_${parcellation_name}_with_19Subcortical.mat
        else
            mv $output_dir/${output_prefix}_lh2lh.mat $output_dir/${output_prefix}_lh2lh_${parcellation_name}_with_19Subcortical.mat
            mv $output_dir/${output_prefix}_lh2rh.mat $output_dir/${output_prefix}_lh2rh_${parcellation_name}_with_19Subcortical.mat
            mv $output_dir/${output_prefix}_rh2rh.mat $output_dir/${output_prefix}_rh2rh_${parcellation_name}_with_19Subcortical.mat
            mv $output_dir/${output_prefix}_lh2subcort.mat $output_dir/${output_prefix}_lh2subcort_${parcellation_name}_with_19Subcortical.mat
            mv $output_dir/${output_prefix}_rh2subcort.mat $output_dir/${output_prefix}_rh2subcort_${parcellation_name}_with_19Subcortical.mat
            mv $output_dir/${output_prefix}_subcort2subcort.mat $output_dir/${output_prefix}_subcort2subcort_${parcellation_name}_with_19Subcortical.mat
            mv $output_dir/${output_prefix}_all2all.mat $output_dir/${output_prefix}_all2all_${parcellation_name}_with_19Subcortical.mat
        endif

        if ( -e $output_dir/${output_prefix}_all2all_${parcellation_name}_with_19Subcortical.mat ) then
            echo "[FC metrics]: Computing ROIs to ROIs correlation matrix sucessed." |& tee -a $LF
        else
            echo "ERROR: Computing ROIs to ROIs correlation matrix failed." |& tee -a $LF
        endif
    endif
endif
echo "=================================== Compute FC metrics finished ==================================" |& tee -a $LF


#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
    echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
    pushd ${CBIG_CODE_DIR}
    git log -1 -- ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/\
CBIG_preproc_FCmetrics_wrapper.csh >> $LF
    popd
endif


exit 0


##################################
##======pass the arguments
##################################
parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
    set flag = $argv[1]; shift;

    switch($flag)
        #subject name
        case "-s":
            if ( $#argv == 0 ) goto arg1err;
            set subject = $argv[1]; shift;
            breaksw
        #path to subject's folder
        case "-d":
            if ( $#argv == 0 ) goto arg1err;
            set sub_dir = $argv[1]; shift;
            breaksw
                        
        case "-bld":
            if ( $#argv == 0 ) goto argerr;
            set bold = ($argv[1]); shift;
            breaksw

        #BOLDbase_suffix
        case "-BOLD_stem":
            if ( $#argv == 0 ) goto arg1err;
            set BOLD_stem = "$argv[1]"; shift;
            breaksw

        case "-SURF_stem":
            if ( $#argv == 0 ) goto arg1err;
            set surf_stem = "$argv[1]"; shift;
            breaksw

        case "-OUTLIER_stem":
            if ( $#argv == 0 ) goto arg1err;
            set outlier_stem = $argv[1]; shift;
            breaksw

        case "-lh_annot":
            if ( $#argv == 0 ) goto arg1err;
            set lh_cortical_ROIs_file = $argv[1]; shift;
            breaksw

        case "-rh_annot":
            if ( $#argv == 0 ) goto arg1err;
            set rh_cortical_ROIs_file = $argv[1]; shift;
            breaksw

        case "-network":
            if ( $#argv == 0 ) goto arg1err;
            set network = "$argv[1]";shift;
            breaksw

        case "-res":
            if ( $#argv == 0 ) goto arg1err;
            set res = "$argv[1]";shift;
            breaksw

        case "-p_type":
            if ( $#argv == 0 ) goto arg1err;
            set parcellation_type = "$argv[1]";shift;
            breaksw

        case "-p_name":
            if ( $#argv == 0 ) goto arg1err;
            set parcellation_name = "$argv[1]";shift;
            breaksw

        case "-censor":
            set censor = 1;
            breaksw

        case "-Pearson_r":
            set Pearson_r = 1;
            breaksw

        case "-nocleanup":
            set nocleanup = 1;
            breaksw

        #running code on MATLAB Runtime (optional)
        case "-matlab_runtime":
            set matlab_runtime = 1;
            breaksw

        default:
            echo ERROR: Flag $flag unrecognized.
            echo $cmdline
            exit 1
            breaksw
    endsw
end
goto parse_args_return;


############################################
##======check passed parameters
############################################
check_params:
if ( "$subject" == "" ) then
    echo "ERROR: subject not specified"
    exit 1;
endif
 
if ( "$sub_dir" == "" ) then
    echo "ERROR: path to subject folder not specified"
    exit 1;
endif

if ( "$bold" == "" ) then
        echo "ERROR: bold number not specified"
        exit 1;
endif

if ( "$BOLD_stem" == "" ) then
    echo "ERROR: BOLD stem not specified"
    exit 1;
endif

if ( "$surf_stem" == "" ) then
    echo "ERROR: SURF stem not specified"
    exit 1;
endif

if ( "$outlier_stem" == "" ) then
    echo "ERROR: outlier stem not specified"
    exit 1;
endif

if ( "$Pearson_r" == 0 ) then
    echo "ERROR: -Pearson_r needed to be passed in." 
    exit 1;
endif

if ( ( "$parcellation_name" == "" && "$lh_cortical_ROIs_file" != "" ) || ( "$parcellation_name" == "" && "$rh_cortical_ROIs_file" != "" ) ) then
    echo "ERROR: parcellation_name not specified"
    exit 1
endif

if ( ( "$parcellation_name" != "" && "$lh_cortical_ROIs_file" == "" ) || ( "$parcellation_name" != "" && "$rh_cortical_ROIs_file" == "" ) ) then
    echo "ERROR: cortical ROIs file not specified"
    exit 1
endif

if ( ( "$network" != "" ) && ( "$network" != "Yeo7" ) && ( "$network" != "Yeo17" ) && ( "$network" != "Kong17" ) ) then
    echo "ERROR: parcellation type not supported, should be 'Yeo7', 'Yeo17' or 'Kong17'"
    exit 1;
endif

if ( ( "$parcellation_type" != "" ) && ( "$parcellation_type" != "Yan" ) && ( "$parcellation_type" != "Schaefer" ) ) then
    echo "ERROR: network type not supported, should be 'Yan' or 'Schaefer'"
    exit 1;
endif

if ( "$res" != "" )  then
    if ( ("$res" > 1000) || ("$res" < 100) || (`expr $res % 100` != 0)) then 
        echo "ERROR: resolution should be from 100 to 1000, multiples of 100"
        exit 1;
    endif
endif

goto check_params_return;


#######################################################
##======Error message
#######################################################
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1

arg2err:
  echo "ERROR: flag $flag requires two arguments"
  exit 1

argerr:
  echo "ERROR: flag $flag requires at least one argument"
  exit 1
  
  
#######################################################
# usage sxit
#######################################################
BEGINHELP

NAME:
    CBIG_preproc_FCmetrics_wrapper.csh

DESCRIPTION:
    Compute FC (functional connectivity) metrics for cortical and 19 subcortical ROIs.
    Parcellation used in this function can be set to Schaefer2018 or Yan2023. Please refer to -p_type, -network, and -res.
    Custom annot files defining cortical ROIs can also be passed in. Please refer to -lh_annot, -rh_annot, and -p_name.

    ***Important***
    We change our default parcellation to Yan2023 homotopic with Kong2022 17 networks, 
    400 parcels resolution. Please refer to -p_type if you want to use Schaefer2018 parcellation.

    The default lh_annot is 
    "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Yan2023_homotopic/
    parcellations/FreeSurfer/fsaverage6/label/kong17/lh.400Parcels_Kong2022_17Networks.annot"
    The default rh_annot is
    "$CBIG_CODE_DIR/CBIG_private/stable_projects/brain_parcellation/Yan2023_homotopic/
    parcellations/FreeSurfer/fsaverage6/label/kong17/rh.400Parcels_Kong2022_17Networks.annot"

    The subcortical ROIs are 19 labels extracted from aseg in subject-specific functional
    space. The 19 labels include:
        8   Left-Cerebellum-Cortex
        10  Left-Thalamus-Proper
        11  Left-Caudate
        12  Left-Putamen
        13  Left-Pallidum
        16  Brain-Stem
        17  Left-Hippocampus
        18  Left-Amygdala
        26  Left-Accumbens-area
        28  Left-VentralDC
        47  Right-Cerebellum-Cortex
        49  Right-Thalamus-Proper
        50  Right-Caudate
        51  Right-Putamen
        52  Right-Pallidum
        53  Right-Hippocampus
        54  Right-Amygdala
        58  Right-Accumbens-area 
        60  Right-VentralDC
    where the numbers are defined by FreeSurfer color table. We choose these 19 subcortical ROIs to be consistent with 
    CIFTI grayordinate format. The final subcortical FC metrics will follow the ascending order of these ROIs 
    (same as this list).

    This function does the following three steps:
        1. Extract the 19 subcortical ROIs from aseg in subject-specific functional space.
        2. Generate input lists for the Matlab function
        3. Call the Matlab function to compute FC metrics for both cortical and subcortical ROIs. 
           Currently we only support compute static Pearson's correlation by using "-Pearson_r" 
           option. In the future we will include other types of metrics.

REQUIRED ARGUMENTS:
    -s  subject :
     fMRI subject id
     
    -d  sub_dir : 
     absolute path to <subject>. All preprocessed data of this subject are assumed to be stored in <sub_dir>/<subject>.
     
    -bld  bold : 
     all bold run numbers of this subject. Each number must has three digits. If this subject has multiple runs, 
     use a space as delimiter, e.g. '002 003'. NOTE: quote sign is necessary.
     
    -BOLD_stem  BOLD_stem :
     stem of input volume. This volume is used to extract the subcortical timeseries in subject-specific functional 
     space. E.g. if the input file name is Sub0001_Ses1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_
     0.009_0.08.nii.gz, then <BOLD_stem> = _rest_stc_mc_resid. This input file is assumed to be stored in <sub_dir>/
     <subject>/bold/<run_number>.
     
    -SURF_stem  surf_stem :
     stem of input surface. The surfaces files are used to extract the cortical timseries. E.g. if the input 
     file names are: lh(/rh).Sub0001_Ses1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_\
     sm6.nii.gz, then <surf_stem> = _rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6. 
     The surface files are assumed to be stored in <sub_dir>/<subject>/surf, and should be in the same space as 
     "lh_cortical_ROIs_file" and "rh_cortical_ROIs_file" (i.e. if your "surf_stem" is xxx_fs5, your 
     "lh_cortical_ROIs_file" and "rh_cortical_ROIs_file" need to be in fsaverage5 space as well).
     
    -OUTLIER_stem  outlier_stem : 
     outlier text file stem. The number of lines in this file is equal to the total number of frames. Each line is 
     a number of 0 or 1 (0-censored, 1-keep). E.g. if the outlier file is Sub0001_Ses1_bld002_FDRMS0.2_DVARS50_motion
     _outliers.txt, then <outlier_stem> = _FDRMS0.2_DVARS50_motion_outliers.txt. The outlier file is assumed to be 
     stored in <sub_dir>/<subject>/qc.
     
    -Pearson_r :
     to compute ROIs to ROIs static Pearson's correlation.
                                  
OPTIONAL ARGUMENTS:
    -p_type parcellation_type
     Parcellation used in computing FC matrix. Currently, we support Schaefer2018 and Yan2023 parcellation.
     Input should be 'Yan' or 'Schaefer' (case sensitive). Default parcellation type is 'Yan'.
     To specify the resolution and matched networks, please refer to -res and -network.
     For custom parcellation, please refer to -lh_annot, -rh_annot and -p_name.

    -res resolution 
     Schaefer2018 and Yan2023 parcellation have different resolutions. Input should be hundreds between 100 to 1000 from 100 to 1000. 
     Default resolution is 400.
     -res only works with Schaefer2018 and Yan2023 parcellation. This argument will be ignored if custom parcellation files are used.

    -network network
     Schaefer2018 and Yan2023 parcellation have different versions, with parcels matched to different versions of networks. 
     Input should be 'Yeo7', 'Yeo17', or 'Kong17' (case sensitive). Default network is 'Kong17'.
     -network only works with Schaefer2018 and Yan2023 parcellation. This argument will be ignored if custom parcellation files are used.

    -lh_annot  lh_cortical_ROIs_file : 
     The annot file including cortical ROIs of the left hemisphere (with absolute path). 
     To use a custom parcellation, -lh_annot, -rh_annot, and -p_name are needed.
     For Schaefer2018 and Yan2023 parcellation, please refer to -p_type, -res, and -network.
     -p_type, -res, and -network will be ignored if -lh_annot, -rh_annot, and -p_name are used.

    -rh_annot  rh_cortical_ROIs_file :
     The annot file including cortical ROIs of the right hemisphere (with absolute path). 
     To use a custom parcellation, -lh_annot, -rh_annot, and -p_name are needed.
     For Schaefer2018 and Yan2023 parcellation, please refer to -p_type, -res, and -network.
     -p_type, -res, and -network will be ignored if -lh_annot, -rh_annot, and -p_name are used.
      
    -p_name parcellation_name
     The parcellation name of your ROIs files. 
     Please note that if parcellation_name is passed in, both lh_cortical_ROIs_file and rh_cortical_ROIs_file should also 
     be passed in.
     Please note that if parcellation_name and cortical_ROIs_file is passed in, -p_type, -network and -res
     will be ignored.
     For Scharfer2018 and Yan2023 cases, please use -p_type, -network and -res, and the parcellation name will be set 
     to <p_type>_<res>Parcels_<network>Networks_order

    -censor :
     ignore censored frames. For example, if the users want to compute ROIs to ROIs correlation (-Pearson_r is used), 
     -censor means the correlations are only computed on uncensored frames.
     
    -nocleanup :
     do not remove intermediate result. For example, if -Pearson_r is used, intermediate files are lh2lh, lh2rh, rh2rh, 
     lh2subcortical, rh2subcortical, and subcortical2subcortical correlation files.

    -matlab_runtime :
     running MATLAB code using MATLAB Runtime instead of MATLAB.

EXAMPLE with default setting:
    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_FCmetrics_wrapper.csh -s Sub0001_Ses1 
    -d ~/storage/fMRI_preprocess -bld "002 003" -BOLD_stem _rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08 
    -SURF_stem _rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6 -OUTLIER_stem 
    _FDRMS0.2_DVARS50_motion_outliers.txt -Pearson_r -censor -nocleanup

EXAMPLE with cortical_ROIs_file and parcellation_name:
    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_FCmetrics_wrapper.csh -s Sub0001_Ses1 
    -d ~/storage/fMRI_preprocess -bld "002 003" -BOLD_stem _rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08 
    -SURF_stem _rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6 -OUTLIER_stem 
    _FDRMS0.2_DVARS50_motion_outliers.txt -Pearson_r -p_name Schaefer2018_400Parcels_17Networks_order 
    -lh_annot 
    "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/
    label/lh.Schaefer2018_400Parcels_17Networks_order.annot" 
    -rh_annot
    "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/
    label/rh.Schaefer2018_400Parcels_17Networks_order.annot"
    -censor -nocleanup

EXAMPLE with different parcellation, resolution and network:
    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_FCmetrics_wrapper.csh -s Sub0001_Ses1 
    -d ~/storage/fMRI_preprocess -bld "002 003" -BOLD_stem _rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08 
    -SURF_stem _rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6 -OUTLIER_stem 
    _FDRMS0.2_DVARS50_motion_outliers.txt -Pearson_r -censor 
    -nocleanup -p_type Schaefer -res 400 -network Yeo17



Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

