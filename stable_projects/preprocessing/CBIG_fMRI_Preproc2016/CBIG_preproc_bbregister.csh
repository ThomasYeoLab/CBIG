#! /bin/csh -f

##########################################
# CBIG preproc bbregister
##########################################
# AUTHOR #################################
# Nanbo Sun 
# 2016/06/18  
##########################################
##########################################
# In this script, we
# 1) Apply bbregister with FSL initialization to each bold run 
# 2) Choose the best run with low FSL cost
#Example: 
#    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_bbregister.csh 
#    -s Sub0001_Ses1 -d ~/storage/fMRI_preprocessing -anat_s Sub0001_Ses1_FS
#    -anat_d ~/storage/sMRI_preprocess -bld '002 003' -BOLD_stem _rest_skip4_stc_mc
##########################################
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


set VERSION = '$Id: CBIG_preproc_bbregister.csh, v 1.0 2016/06/18 $'

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
set subject_dir = ""
set anat = ""
set zpdbold = ""
set BOLD_stem = ""
set force = 0        # Default if file exist, skip the step

set root_dir = `python -c "import os; print(os.path.realpath('$0'))"`
set root_dir = `dirname $root_dir`

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

cd $subject_dir/$subject

#############################################
# BOLD Information
#############################################

set qc = $subject_dir/$subject"/qc"

if (! -e $qc) then
     mkdir -p $qc
endif
if (! -e logs) then
     mkdir -p logs
endif
set LF = $subject_dir/$subject/logs/CBIG_preproc_bbregister.log
if( -e $LF ) then
    rm $LF
endif
touch $LF
echo "[REG]: logfile = $LF"
echo "Registration" >> $LF
echo "[CMD]: CBIG_preproc_bbregister.csh $cmdline"   >>$LF

set boldfolder = "$subject_dir/$subject/bold"
echo "[REG]: boldfolder = $boldfolder" |& tee -a $LF
echo "[REG]: zpdbold = $zpdbold" |& tee -a $LF


#############################################
# Apply bbregister with FSL initialization to each bold run
#############################################

echo "=======================bbregister with FSL initialization=======================" |& tee -a $LF
cd $boldfolder

foreach curr_bold ($zpdbold)
    pushd $curr_bold
    mkdir -p bbr_orig
    set boldfile = $subject"_bld"$curr_bold$BOLD_stem
    if ( (! -e bbr_orig/$boldfile'_reg.dat') || ($force == 1)) then
        echo "[REG]: boldfile = $boldfile" |& tee -a $LF
        set cmd = "bbregister --bold --s $anat --init-fsl --mov $boldfile.nii.gz --reg bbr_orig/$boldfile'_reg.dat'"
        echo $cmd |& tee -a $LF
        eval $cmd >> $LF
    else
        echo "bbr_orig/$boldfile'_reg.dat' already exists"|& tee -a $LF 
    endif
    rsync -az bbr_orig/$boldfile"_reg.dat" $boldfile"_reg.dat" 
    popd
end
echo "=======================FSL initialization done!=======================" |& tee -a $LF
echo "" |& tee -a $LF


#############################################
# Choose the best run with low FSL cost
#############################################

echo "=======================choose the best run =======================" |& tee -a $LF
# grab registration cost function values
set reg_cost_file = $qc/CBIG_preproc_bbregister_intra_sub_reg.cost
if (-e $reg_cost_file) then
    rm $reg_cost_file
endif
foreach curr_bold ($zpdbold)
    set boldfile = $subject"_bld"$curr_bold$BOLD_stem
    cat $curr_bold/bbr_orig/$boldfile'_reg.dat.mincost' | awk '{print $1}' >> $reg_cost_file
end

# compute best fsl cost
set init_fsl = `cat $reg_cost_file`
set min_fsl_cost = 100000
set count = 1;
while ($count <= $#init_fsl)
    set comp = `echo "$init_fsl[$count] < $min_fsl_cost" | bc`
    if ($comp == 1) then
        set best_fsl_index = $count
        set min_fsl_cost = $init_fsl[$count]
    endif
    @ count = $count + 1;
end
echo "Best fsl register is run $zpdbold[$best_fsl_index] with cost = $min_fsl_cost" |& tee -a $LF
set best_run = $zpdbold[$best_fsl_index]

# save best run number in the qc folder
set best_run_file = $qc/CBIG_preproc_bbregister_best_run.dat
if (-e $best_run_file) then
    rm $best_run_file
endif
echo $best_run >> $best_run_file

# check if the best registration reduce the bbr cost of other run. If so, use that registration
set bestboldfile = $subject"_bld"$best_run$BOLD_stem
foreach curr_bold ($zpdbold)
    if ($curr_bold != $best_run) then
        set boldfile = $subject"_bld"$curr_bold$BOLD_stem
        set cmd = "bbregister --bold --s $anat --init-reg $best_run/bbr_orig/$bestboldfile'_reg.dat'"
        set cmd = "$cmd --mov $curr_bold/$boldfile.nii.gz --reg $curr_bold/bbr_use_best_run/$boldfile'_reg.dat'"
        set cmd = "$cmd --initcost $curr_bold/bbr_use_best_run/$boldfile'_use_best_run_initcost.dat'"
        echo $cmd |& tee -a $LF
        eval $cmd
        set use_best_run_cost = `cat $curr_bold/bbr_use_best_run/$boldfile'_use_best_run_initcost.dat' | awk '{print $1}'`
        set original_cost = `cat $curr_bold/bbr_orig/$boldfile'_reg.dat.mincost' | awk '{print $1}'`
        echo "[REG]: BBR cost of run $curr_bold using transformation matrix from the best run: $use_best_run_cost" >> $LF
        echo "[REG]: BBR cost of run $curr_bold using its own transformation matrix: $original_cost" |& tee -a $LF
        set comp = `echo "$use_best_run_cost < $original_cost" | bc`
        if ($comp == 1) then
            echo "[Reg] Registration from best run (run$best_run) reduces the bbr cost of run$curr_bold" |& tee -a $LF
            echo "[Reg]: Registration from the best run will be applied to run$curr_bold" |& tee -a $LF
            set cmd = "rsync -az $best_run/bbr_orig/$bestboldfile'_reg.dat' $curr_bold/$boldfile'_reg.dat'"
            echo $cmd |& tee -a $LF
            eval $cmd
            set cmd = "rsync -az $curr_bold/bbr_use_best_run/$boldfile'_use_best_run_initcost.dat'"
            set cmd = "$cmd $curr_bold/$boldfile'_reg.dat.mincost'"
            echo $cmd |& tee -a $LF
            eval $cmd
        else
            echo "[Reg]: Registration from best run (run$best_run) doesn't reduces the bbr cost of run$curr_bold" |& tee -a $LF
            echo "[REG]: BBR cost of run $curr_bold will keep using its original transformation matrix" |& tee -a $LF
            set cmd = "rsync -az $curr_bold/bbr_orig/$boldfile'_reg.dat.mincost' $curr_bold/$boldfile'_reg.dat.mincost'"
            echo $cmd |& tee -a $LF
            eval $cmd
        endif
        if ($boldfolder != "" && $curr_bold != "") then
            rm $boldfolder/$curr_bold/bbr_use_best_run/${subject}*
            rmdir $boldfolder/$curr_bold/bbr_use_best_run/
        endif
    endif
end
#copy the bbr cost of the best run to its bold run folder
set cmd = "rsync -az $best_run/bbr_orig/$bestboldfile'_reg.dat.mincost' $best_run/$bestboldfile'_reg.dat.mincost'"
echo $cmd |& tee -a $LF
eval $cmd

#update bbr cost in the qc folder
rm $reg_cost_file
foreach curr_bold ($zpdbold)
    set boldfile = $subject"_bld"$curr_bold$BOLD_stem
    cat $curr_bold/$boldfile'_reg.dat.mincost' | awk '{print $1}' >> $reg_cost_file
end

#########################
# Create loose brain mask and apply it to current fMRI volumes
#########################
echo "=======================Create loose whole brain mask for saving space=======================" |& tee -a $LF
set REG_stem = $BOLD_stem"_reg"
set MASK_stem = $BOLD_stem
echo $root_dir |& tee -a $LF
set cmd = "$root_dir/CBIG_preproc_create_mask.csh -s $subject -d $subject_dir -anat_s $anat -anat_d $anat_dir"
set cmd = "$cmd -bld '$zpdbold' -REG_stem $REG_stem -MASK_stem $MASK_stem -loose_whole_brain"
echo $cmd |& tee -a $LF
eval $cmd

echo "=======================Apply loose whole brain mask to current fMRI volumes=======================" |& tee -a $LF
foreach curr_bold ($zpdbold)
    set boldfile = $subject"_bld"$curr_bold$BOLD_stem
    set cmd = "fslmaths $curr_bold/$boldfile -mas mask/$subject.loosebrainmask.bin.nii.gz $curr_bold/$boldfile"
    echo $cmd |& tee -a $LF
    eval $cmd
end

#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
    echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
    pushd ${CBIG_CODE_DIR}
    git log -1 -- ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_bbregister.csh \
    >> $LF
    popd
endif

echo "****************************************************************" |& tee -a $LF
exit 0;

##########################################
# Parse Arguments 
##########################################

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
            set subject_dir = $argv[1]; shift;
            breaksw

        #anatomical ID
        case "-anat_s":
            if ($#argv == 0) goto arg1err;
            set anat = $argv[1]; shift;
            breaksw

        #anatomical directory
        case "-anat_d":
            if ($#argv == 0) goto arg1err;
            set anat_dir = $argv[1];
            setenv SUBJECTS_DIR $argv[1]; shift;
            breaksw

        #bold run number
        case "-bld":
            if ( $#argv == 0 ) goto arg1err;
            set zpdbold = ($argv[1]); shift;
            breaksw

        #input file stem
        case "-BOLD_stem":
            if ( $#argv == 0 ) goto arg1err;
            set BOLD_stem = $argv[1]; shift;
            breaksw

        #update results, if exist then overwrite
        case "-force":
            set force = 1;
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

if ( "$subject" == "" ) then
    echo "ERROR: subject not specified"
    exit 1;
endif
if ( "$subject_dir" == "" ) then
    echo "ERROR: path to subject folder not specified"
    exit 1;
endif
if ( "$anat" == "" ) then
    echo "ERROR: anatomical folder not specified"
    exit 1;
endif
if ( "$zpdbold" == "" ) then
    echo "ERROR: bold run not specified"
    exit 1;
endif
if ( "$BOLD_stem" == "" ) then
    echo "ERROR: input file stem not specified"
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
BEGINHELP

NAME:
    CBIG_preproc_bbregister.csh

DESCRIPTION:
    T1-T2* registration using bbregister.
    We use the bbregister default parameters. For FreeSurfer 5.3, the reference frame is 0, 
    i.e. the 1st frame.

    This function 
    1) Applies bbregister with FSL initialization to each bold run. 
    2) Chooses the best run with lowest bbr cost
    3) Applies the registration matrix from best run to all other runs and computes the bbr cost.
    Therefore, for each run, we have two bbr costs: one cost computed using original regitration 
    matrix from step 1 and one cost computed using the best run registration matrix.
    Then compare the two cost. If the cost using best run registration is lower, the final
    registration of this run 
    <subject_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_reg.dat
    will be the best run registration and the final cost
    <subject_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_reg.dat.mincost
    will be the cost computed using the best run registration. Else, they will be the original 
    registration and cost from step 1
    4) Create a loose whole brain mask using the bold run with lowest bbr cost and apply the maske to 
    fMRI volumes that were used to do bbregister. 
    PS: The aim of Step4 is to save space because .nii.gz file will compress more if many voxels in the
    volume are 0s.

REQUIRED ARGUMENTS:
    -s  subject_id           : name of the subject
    -d  subject_dir          : absolute path to <subject_id>, i.e. all preprocessed data 
                               are assumed to be stored in <subject_dir>/<subject_id>.
    -anat_s  anat_src        : name of anatomical folder for the subject, i.e. recon-all results.
    -anat_d  anat_dir        : absolute path to <anat_dir>, i.e. all recon-all results are 
                               assumed to be stored in <anat_dir>/<anat_src>.
    -bld  bold_runs          : all bold run numbers of this subject. Each number must be three 
                               digits. If this subject has multiple runs, please use space as 
                               delimiter between two runs (e.g. -bld '002 003'). NOTE: quote sign 
                               is necessary.
    -BOLD_stem  BOLD_stem    : stem of input file. E.g. if the input file is 
                               Sub0001_Ses1_bld002_rest_skip4_stc_mc.nii.gz, then <BOLD_stem> = 
                               _rest_skip4_stc_mc. This input file is assumed to be in 
                               <subject_dir>/<subject_id>/bold/<run_number>.

OPTIONAL ARGUMENTS:
    -force                   : update results, if exist then overwrite
    -help                    : help
    -version                 : version

OUTPUTS:
    Files generated by bbregister:
    <subject_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_reg.dat
    <subject_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_reg.dat.mincost

EXAMPLE:
    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_bbregister.csh 
    -s Sub0001_Ses1 -d ~/storage/fMRI_preprocessing -anat_s Sub0001_Ses1_FS
    -anat_d ~/storage/sMRI_preprocess -bld '002 003' -BOLD_stem _rest_skip4_stc_mc
