#! /bin/csh
#

####################################################
# export paths for subdirectories of CBIG repository
####################################################
setenv CBIG_BIN_DIR $CBIG_CODE_DIR/bin
setenv CBIG_SD_DIR  $CBIG_CODE_DIR/external_packages/SD/SDv1.5.1-svn593
set path = ($CBIG_BIN_DIR $path)

##################
# setup FreeSurfer
##################

if ((! $?FREESURFER_HOME)  || { eval 'test ! -d $FREESURFER_HOME' }) then
  if ($?prompt) then
    echo "WARNING: FREESURFER_HOME is not set or points to a non-existing directory"
  endif
else
  setenv FSF_OUTPUT_FORMAT nii.gz
  setenv FSFAST_HOME $FREESURFER_HOME/fsfast
  
  if ($?prompt) then
    source $FREESURFER_HOME/SetUpFreeSurfer.csh
  else
    source $FREESURFER_HOME/SetUpFreeSurfer.csh > /dev/null
  endif

  # check whether user's default freesurfer version matches CBIG default freesurfer version
  set CBIG_default_FREESURFER = `cat $CBIG_CODE_DIR/setup/CBIG_sample_config.csh | grep FREESURFER_HOME \
    | grep -E -o "[0-9]+.[0-9]+.[0-9]+"`
  set USER_default_FREESURFER = `freesurfer | grep -E "v[0-9]+.[0-9]+.[0-9]+" | grep -E -o "[0-9]+.[0-9]+.[0-9]+"`
  
  if (! $?USER_default_FREESURFER) then
    if ($?prompt) then
      echo "[WARNING]: This version of CBIG repository's default FREESURFER version is $CBIG_default_FREESURFER."
      echo -n "You appear to be using an unrecognized/invalid FREESURFER version."
      echo "Switch to FREESURFER version $CBIG_default_FREESURFER if possible.\n"
    endif
  else if ("$USER_default_FREESURFER" != "$CBIG_default_FREESURFER") then
    if ($?prompt) then
      echo "[WARNING]: This version of CBIG repository's default FREESURFER version is $CBIG_default_FREESURFER."
      echo -n "You appear to be using FREESURFER $USER_default_FREESURFER (which may or may not work with current repo)."
      echo "Switch to FREESURFER version $CBIG_default_FREESURFER if possible."
      echo -n "Note: stable projects may not follow default setting, "
      echo "refer to the proper config file of the project you use."
      echo -n "If you are unsure, please run our stable projects' examples "
      echo "to check if the results can be replicated.\n"
    endif
  endif
endif

###############
# set up Matlab
###############

if ((! $?CBIG_MATLAB_DIR)  || { eval 'test ! -d $CBIG_MATLAB_DIR' }) then
  if ($?prompt) then
    echo "WARNING: CBIG_MATLAB_DIR is not set or points to a non-existing directory"
  endif
else
  set path = ($CBIG_MATLAB_DIR/bin $path)

  # check whether user's default matlab version matches CBIG default matlab version
  set CBIG_default_MATLAB = `cat $CBIG_CODE_DIR/setup/CBIG_sample_config.csh | grep -o "R20[0-9][0-9][a-b]"`
  set USER_default_MATLAB = `which matlab | sed 's/bin\/matlab//g'`
  set USER_default_MATLAB = `cat $USER_default_MATLAB/help/matlab/matlab-version-and-license.html \
    | grep -o "R20[0-9][0-9][a-b]" | sed -n '1p'`

  if (! $?USER_default_MATLAB) then
    if ($?prompt) then
      echo "[WARNING]: This version of CBIG repository's default MATLAB version is $CBIG_default_MATLAB."
      echo -n "You appear to be using an unrecognized/invalid MATLAB version."
      echo "Switch to MATLAB $CBIG_default_MATLAB if possible.\n"
    endif
  else if ("$USER_default_MATLAB" != "$CBIG_default_MATLAB") then
    
      echo "[WARNING]: This version of CBIG repository's default MATLAB version is $CBIG_default_MATLAB."
      echo -n "You appear to be using MATLAB $USER_default_MATLAB (which may or may not work with current repo)."
      echo "Switch to MATLAB version $CBIG_default_MATLAB if possible."
      echo -n "Note: stable projects may not follow default setting, "
      echo "refer to the proper config file of the project you use."
      echo -n "If you are unsure, please run our stable projects' examples "
      echo "to check if the results can be replicated.\n"
    
  endif
endif

############
# set up FSL
############

if ((! $?CBIG_FSLDIR)  || { eval 'test ! -d $CBIG_FSLDIR' }) then
  if ($?prompt) then
    echo "WARNING: CBIG_FSLDIR is not set or points to a non-existing directory"
  endif
else
  setenv FSLDIR $CBIG_FSLDIR #used by FSL
  setenv FSL_DIR $FSLDIR #used by Freesurfer
  source $FSLDIR/etc/fslconf/fsl.csh
  set path = ($FSL_DIR/bin $path)

  if ($FSL_DIR != $FSLDIR) then
    if ($?prompt) then
      echo "WARNING: FSL_DIR = $FSL_DIR"
      echo "WARNING: FSLDIR = $FSLDIR"
      echo "WARNING: FSLDIR and FSL_DIR doesn't match! This may cause some problem! \
FSLDIR is the FSL directory defined by FSL. FSL_DIR is the FSL directory defined by FreeSurfer."
    endif
  endif

  # check whether user's default fsl version matches CBIG default fsl version
  set CBIG_default_FSL = `cat $CBIG_CODE_DIR/setup/CBIG_sample_config.csh | grep CBIG_FSLDIR \
    | grep -E -o "[0-9]+.[0-9]+.[0-9]+"`
  set USER_default_FSL = `which fsl | sed 's/bin\/fsl//g'`
  set USER_default_FSL = `cat $USER_default_FSL/etc/fslversion`

  if (! $?USER_default_FSL) then
    if ($?prompt) then
      echo "[WARNING]: This version of CBIG repository's default FSL version is $CBIG_default_FSL."
      echo -n "You appear to be using an unrecognized/invalid FSL version."
      echo "Switch to FSL $CBIG_default_FSL if possible.\n"
    endif
  else if ("$USER_default_FSL" != "$CBIG_default_FSL") then
    if ($?prompt) then
      echo "[WARNING]: This version of CBIG repository's default FSL version is $CBIG_default_FSL."
      echo -n "You appear to be using FSL $USER_default_FSL (which may or may not work with current repo)."
      echo "Switch to FSL version $CBIG_default_FSL if possible."
      echo -n "Note: stable projects may not follow default setting, "
      echo "refer to the proper config file of the project you use."
      echo -n "If you are unsure, please run our stable projects' examples "
      echo "to check if the results can be replicated.\n"
    endif
  endif
endif

#############
# set up ANTS
#############

if ((! $?CBIG_ANTS_DIR)  || { eval 'test ! -d $CBIG_ANTS_DIR' }) then
  if ($?prompt) then
    echo "WARNING: CBIG_ANTS_DIR is not set or points to a non-existing directory"
  endif
else
  setenv ANTSPATH $CBIG_ANTS_DIR
  set path = ($CBIG_ANTS_DIR $path)
endif

############################
# set up Workbench directory
############################

if ((! $?CBIG_WB_DIR)  || { eval 'test ! -d $CBIG_WB_DIR' }) then
  if ($?prompt) then
    echo "WARNING: CBIG_WB_DIR is not set or points to a non-existing directory"
  endif
else
  set path = ($CBIG_WB_DIR/bin_rh_linux64 $path)

  # check whether user's default workbench version matches CBIG default workbench version
  set CBIG_default_WB = `cat $CBIG_CODE_DIR/setup/CBIG_sample_config.csh | grep CBIG_WB_DIR \
    | grep -E -o "[0-9]+.[0-9]+.[0-9]+"`
  set USER_default_WB = `wb_command | grep "Version" | grep -E -o "[0-9]+.[0-9]+.[0-9]+"`

  if (! $?USER_default_WB) then
    if ($?prompt) then
      echo "[WARNING]: This version of CBIG repository's default WORKBENCH version is $CBIG_default_WB."
      echo -n "You appear to be using an unrecognized/invalid WORKBENCH version."
      echo "Switch to WORKBENCH $CBIG_default_WB if possible.\n"
    endif
  else if ("$USER_default_WB" != "$CBIG_default_WB") then
    if ($?prompt) then
      echo "[WARNING]: This version of CBIG repository's default WORKBENCH version is $CBIG_default_WB."
      echo -n "You appear to be using WORKBENCH $USER_default_WB (which may or may not work with current repo)."
      echo "Switch to WORKBENCH version $CBIG_default_WB if possible."
      echo -n "Note: stable projects may not follow default setting, "
      echo "refer to the proper config file of the project you use."
      echo -n "If you are unsure, please run our stable projects' examples "
      echo "to check if the results can be replicated.\n"
    endif
  endif
endif

#############
# set up AFNI
#############

if ((! $?CBIG_AFNI_DIR)  || { eval 'test ! -d $CBIG_AFNI_DIR' }) then
  if ($?prompt) then
    echo "WARNING: CBIG_AFNI_DIR is not set or points to a non-existing directory"
  endif
else
  set path = ($CBIG_AFNI_DIR $path)
  
  # check whether user's default afni version matches CBIG default afni version
  set CBIG_default_AFNI = `cat $CBIG_CODE_DIR/setup/CBIG_sample_config.csh | grep CBIG_AFNI_DIR \
    | grep -o "afni.*" | grep -o "AFNI.*" | sed 's/\/linux_openmp_64//g'`
  set USER_default_AFNI = `afni --version | grep -o "AFNI_.*" | sed 's/)//g'`

  if (! $?USER_default_AFNI) then
    if ($?prompt) then
      echo "[WARNING]: This version of CBIG repository's default AFNI version is $CBIG_default_AFNI."
      echo -n "You appear to be using an unrecognized/invalid AFNI version."
      echo "Switch to AFNI $CBIG_default_AFNI if possible. \n"
    endif
  else if ("$USER_default_AFNI" != "$CBIG_default_AFNI") then
    if ($?prompt) then
      echo "[WARNING]: This version of CBIG repository's default AFNI version is $CBIG_default_AFNI."
      echo -n "You appear to be using AFNI $USER_default_AFNI (which may or may not work with current repo)."
      echo "Switch to AFNI version $CBIG_default_AFNI if possible."
      echo -n "Note: stable projects may not follow default setting, "
      echo "refer to the proper config file of the project you use. "
      echo -n "If you are unsure, please run our stable projects' examples "
      echo "to check if the results can be replicated.\n"
    endif
  endif
endif

####################################
# create symlink to git hook scripts
####################################

if (! -f "$CBIG_CODE_DIR/.git/hooks/pre-commit") then
  mkdir -p $CBIG_CODE_DIR/.git/hooks
  ln -s "$CBIG_CODE_DIR/hooks/pre-commit" "$CBIG_CODE_DIR/.git/hooks/pre-commit"
endif
if (! -f "$CBIG_CODE_DIR/.git/hooks/pre-push") then
  mkdir -p $CBIG_CODE_DIR/.git/hooks
  ln -s "$CBIG_CODE_DIR/hooks/pre-push" "$CBIG_CODE_DIR/.git/hooks/pre-push"
endif
