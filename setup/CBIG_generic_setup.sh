#!/bin/bash
#

# export paths for subdirectories of CBIG repository
export CBIG_BIN_DIR=$CBIG_CODE_DIR/bin
export CBIG_SD_DIR=$CBIG_CODE_DIR/external_packages/SD/SDv1.5.1-svn593
export PATH=$CBIG_BIN_DIR:$PATH

# set up Matlab
if [ -z $CBIG_MATLAB_DIR ] || [ ! -d $CBIG_MATLAB_DIR ]; then
  echo "[WARNING]: CBIG_MATLAB_DIR is not set or points to a non-existent directory"
else
  export PATH=$CBIG_MATLAB_DIR/bin:$PATH
fi

# set up FSL
if [ -z $CBIG_FSLDIR ] || [ ! -d $CBIG_FSLDIR ]; then
  echo "[WARNING]: CBIG_FSLDIR is not set or points to a non-existent directory"
else
  export FSLDIR=$CBIG_FSLDIR #used by FSL
  export FSL_DIR=$FSLDIR #used by Freesurfer
  source $FSLDIR/etc/fslconf/fsl.sh
  export PATH=$FSL_DIR/bin:$PATH
  
  if [ "$FSL_DIR" != "$FSLDIR" ]
  then
    echo "[WARNING]: FSL_DIR = $FSL_DIR"
    echo "[WARNING]: FSLDIR = $FSLDIR"
    echo "[WARNING]: FSLDIR and FSL_DIR doesn't match! This may cause some problem! FSLDIR is the FSL directory defined by FSL. FSL_DIR is the FSL directory defined by FreeSurfer."
  fi
fi

# setup FreeSurfer
if [ -z $FREESURFER_HOME ] || [ ! -d $FREESURFER_HOME ]; then
  echo "[WARNING]: FREESURFER_HOME is not set or points to a non-existent directory"
else
  export FSF_OUTPUT_FORMAT="nii.gz"
  export FSFAST_HOME=$FREESURFER_HOME/fsfast
  dumb_variable="dumb"
  if [ "${term}" == "${dumb_variable}" ]
  then
    source $FREESURFER_HOME/SetUpFreeSurfer.sh > /dev/null
  else
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
  fi
fi

# set up ANTS
if [ -z $CBIG_ANTS_DIR ] || [ ! -d $CBIG_ANTS_DIR ]; then
  echo "[WARNING]: CBIG_ANTS_DIR is not set or points to a non-existent directory"
else
  export ANTSPATH=$CBIG_ANTS_DIR
  export PATH=$CBIG_ANTS_DIR:$PATH
fi

# set up Workbench directory
if [ -z $CBIG_WB_DIR ] || [ ! -d $CBIG_WB_DIR ]; then
  echo "[WARNING]: CBIG_WB_DIR is not set or points to a non-existent directory"
else
  export PATH=$CBIG_WB_DIR/bin_rh_linux64:$PATH
fi

# setup Caret directory
if [ -z $CBIG_CARET_DIR ] || [ ! -d $CBIG_CARET_DIR ]; then
  echo "[WARNING]: CBIG_CARET_DIR is not set or points to a non-existent directory"
else
  export PATH=$CBIG_CARET_DIR/bin_linux64:$PATH
fi

# set up AFNI
if [ -z $CBIG_AFNI_DIR ] || [ ! -d $CBIG_AFNI_DIR ]; then
  echo "[WARNING]: CBIG_AFNI_DIR is not set or points to a non-existent directory"
else
  export PATH=$CBIG_AFNI_DIR:$PATH
fi



# create symlink to git hook scripts
if [ ! -f "$CBIG_CODE_DIR/.git/hooks/pre-commit" ];
then
  ln -s "$CBIG_CODE_DIR/hooks/pre-commit" "$CBIG_CODE_DIR/.git/hooks/pre-commit"
fi
if [ ! -f "$CBIG_CODE_DIR/.git/hooks/pre-push" ];
then
  ln -s "$CBIG_CODE_DIR/hooks/pre-push" "$CBIG_CODE_DIR/.git/hooks/pre-push"
fi
