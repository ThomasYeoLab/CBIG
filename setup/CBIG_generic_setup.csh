#! /bin/csh
#

# export paths for subdirectories of CBIG repository
setenv CBIG_BIN_DIR $CBIG_CODE_DIR/bin
setenv CBIG_SD_DIR  $CBIG_CODE_DIR/external_packages/SD/SDv1.5.1-svn593
set path = ($CBIG_BIN_DIR $path)

# set up Matlab
if ((! $?CBIG_MATLAB_DIR)  || { eval 'test ! -d $CBIG_MATLAB_DIR' }) then
  echo "WARNING: CBIG_MATLAB_DIR is not set or points to a non-existing directory"
else
  set path = ($CBIG_MATLAB_DIR/bin $path)
endif

# set up FSL
if ((! $?CBIG_FSLDIR)  || { eval 'test ! -d $CBIG_FSLDIR' }) then
  echo "WARNING: CBIG_FSLDIR is not set or points to a non-existing directory"
else
  setenv FSLDIR $CBIG_FSLDIR #used by FSL
  setenv FSL_DIR $FSLDIR #used by Freesurfer
  source $FSLDIR/etc/fslconf/fsl.csh
  set path = ($FSL_DIR/bin $path)

  if ($FSL_DIR != $FSLDIR) then
    echo "WARNING: FSL_DIR = $FSL_DIR"
    echo "WARNING: FSLDIR = $FSLDIR"
    echo "WARNING: FSLDIR and FSL_DIR doesn't match! This may cause some problem! FSLDIR is the FSL directory defined by FSL. FSL_DIR is the FSL directory defined by FreeSurfer."
  endif
endif

# FreeSurfer
if ((! $?FREESURFER_HOME)  || { eval 'test ! -d $FREESURFER_HOME' }) then
  echo "WARNING: FREESURFER_HOME is not set or points to a non-existing directory"
else
  setenv FSF_OUTPUT_FORMAT nii.gz
  setenv FSFAST_HOME $FREESURFER_HOME/fsfast
  if("$term" == "dumb") then
    source $FREESURFER_HOME/SetUpFreeSurfer.csh > /dev/null
  else
    source $FREESURFER_HOME/SetUpFreeSurfer.csh
  endif
endif

# set up ANTS
if ((! $?CBIG_ANTS_DIR)  || { eval 'test ! -d $CBIG_ANTS_DIR' }) then
  echo "WARNING: CBIG_ANTS_DIR is not set or points to a non-existing directory"
else
  setenv ANTSPATH $CBIG_ANTS_DIR
  set path = ($CBIG_ANTS_DIR $path)
endif

# set up Workbench
if ((! $?CBIG_WB_DIR)  || { eval 'test ! -d $CBIG_WB_DIR' }) then
  echo "WARNING: CBIG_WB_DIR is not set or points to a non-existing directory"
else
  set path = ($CBIG_WB_DIR/bin_rh_linux64 $path)
endif

# set up Caret
if ((! $?CBIG_CARET_DIR)  || { eval 'test ! -d $CBIG_CARET_DIR' }) then
  echo "WARNING: CBIG_CARET_DIR is not set or points to a non-existing directory"
else
  set path = ($CBIG_CARET_DIR/bin_linux64 $path)
endif

# set up AFNI
if ((! $?CBIG_AFNI_DIR)  || { eval 'test ! -d $CBIG_AFNI_DIR' }) then
  echo "WARNING: CBIG_AFNI_DIR is not set or points to a non-existing directory"
else
  set path = ($CBIG_AFNI_DIR $path)
endif



# create symlink to git hook scripts
if (! -f "$CBIG_CODE_DIR/.git/hooks/pre-commit") then
  ln -s "$CBIG_CODE_DIR/hooks/pre-commit" "$CBIG_CODE_DIR/.git/hooks/pre-commit"
endif
if (! -f "$CBIG_CODE_DIR/.git/hooks/pre-push") then
  ln -s "$CBIG_CODE_DIR/hooks/pre-push" "$CBIG_CODE_DIR/.git/hooks/pre-push"
endif
