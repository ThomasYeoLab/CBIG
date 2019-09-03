#!/bin/bash
#

##########################
# define a helper function
##########################

warn_msg () {
  # only output warning msg when the shell is interactive
  if [[ $- == *i* ]]; then
    echo $1
  fi
}

####################################################
# export paths for subdirectories of CBIG repository
####################################################

export CBIG_BIN_DIR=$CBIG_CODE_DIR/bin
export CBIG_SD_DIR=$CBIG_CODE_DIR/external_packages/SD/SDv1.5.1-svn593
export PATH=$CBIG_BIN_DIR:$PATH

##################
# setup FreeSurfer
##################

if [ -z $FREESURFER_HOME ] || [ ! -d $FREESURFER_HOME ]; then
  warn_msg
  warn_msg "[WARNING]: FREESURFER_HOME is not set or points to a non-existent directory"
else
  export FSF_OUTPUT_FORMAT="nii.gz"
  export FSFAST_HOME=$FREESURFER_HOME/fsfast

  if [[ $- == *i* ]]; then
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
  else
    source $FREESURFER_HOME/SetUpFreeSurfer.sh > /dev/null
  fi

  # check whether user's default freesurfer version matches CBIG default freesurfer version
  CBIG_default_FREESURFER=`cat $CBIG_CODE_DIR/setup/CBIG_sample_config.sh | grep FREESURFER_HOME \
    | grep -E -o [0-9]+.[0-9]+.[0-9]+`
  USER_default_FREESURFER=`freesurfer | grep -E v[0-9]+.[0-9]+.[0-9]+ | grep -E -o [0-9]+.[0-9]+.[0-9]+`

  if [ -z "$USER_default_FREESURFER" ]; then
    warn_msg
    warn_msg "[WARNING]: This version of CBIG repository's default FREESURFER version is $CBIG_default_FREESURFER."
    warn_msg "You appear to be using an unrecognized/invalid FREESURFER version. \
Switch to FREESURFER version $CBIG_default_FREESURFER if possible."
  elif [ "$USER_default_FREESURFER" != "$CBIG_default_FREESURFER" ]; then
    warn_msg
    warn_msg "[WARNING]: This version of CBIG repository's default FREESURFER version is $CBIG_default_FREESURFER."
    warn_msg "You appear to be using FREESURFER $USER_default_FREESURFER (which may or may not work with current repo).\
 Switch to FREESURFER version $CBIG_default_FREESURFER if possible."
    warn_msg "Note: stable projects may not follow default setting, \
refer to the proper config file of the project you use."
    warn_msg "If you are unsure, please run our stable projects' examples \
to check if the results can be replicated."
  fi
fi

###############
# set up Matlab
###############

if [ -z $CBIG_MATLAB_DIR ] || [ ! -d $CBIG_MATLAB_DIR ]; then
  warn_msg
  warn_msg "[WARNING]: CBIG_MATLAB_DIR is not set or points to a non-existent directory"
else
  export PATH=$CBIG_MATLAB_DIR/bin:$PATH

  # check whether user's default matlab version matches CBIG default matlab version
  CBIG_default_MATLAB=`cat $CBIG_CODE_DIR/setup/CBIG_sample_config.sh | grep -o R20[0-9][0-9][a-b]`
  USER_default_MATLAB=`which matlab | sed 's/bin\/matlab//g'`
  USER_default_MATLAB=`cat $USER_default_MATLAB/help/matlab/matlab-version-and-license.html \
    | grep -o R20[0-9][0-9][a-b] | sed -n '1p'`

  if [ -z "$USER_default_MATLAB" ]; then
    warn_msg
    warn_msg "[WARNING]: This version of CBIG repository's default MATLAB version is $CBIG_default_MATLAB."
    warn_msg "You appear to be using an unrecognized/invalid MATLAB version. \
Switch to MATLAB $CBIG_default_MATLAB if possible."
  elif [ "$USER_default_MATLAB" != "$CBIG_default_MATLAB" ]; then
    warn_msg
    warn_msg "[WARNING]: This version of CBIG repository's default MATLAB version is $CBIG_default_MATLAB."
    warn_msg "You appear to be using MATLAB $USER_default_MATLAB (which may or may not work with current repo). \
Switch to MATLAB version $CBIG_default_MATLAB if possible."
    warn_msg "Note: stable projects may not follow default setting, \
refer to the proper config file of the project you use."
    warn_msg "If you are unsure, please run our stable projects' examples \
to check if the results can be replicated."
  fi
fi

############
# set up FSL
############

if [ -z $CBIG_FSLDIR ] || [ ! -d $CBIG_FSLDIR ]; then
  warn_msg
  warn_msg "[WARNING]: CBIG_FSLDIR is not set or points to a non-existent directory"
else
  export FSLDIR=$CBIG_FSLDIR #used by FSL
  export FSL_DIR=$FSLDIR #used by Freesurfer
  source $FSLDIR/etc/fslconf/fsl.sh
  export PATH=$FSL_DIR/bin:$PATH
  
  if [ "$FSL_DIR" != "$FSLDIR" ]; then
    warn_msg
    warn_msg "[WARNING]: FSL_DIR = $FSL_DIR"
    warn_msg "[WARNING]: FSLDIR = $FSLDIR"
    warn_msg "[WARNING]: FSLDIR and FSL_DIR doesn't match! This may cause some problem! \
      FSLDIR is the FSL directory defined by FSL. FSL_DIR is the FSL directory defined by FreeSurfer."
  fi

  # check whether user's default fsl version matches CBIG default fsl version
  CBIG_default_FSL=`cat $CBIG_CODE_DIR/setup/CBIG_sample_config.sh | grep CBIG_FSLDIR \
    | grep -E -o [0-9]+.[0-9]+.[0-9]+`
  USER_default_FSL=`which fsl | sed 's/bin\/fsl//g'`
  USER_default_FSL=`cat $USER_default_FSL/etc/fslversion`

  if [ -z "$USER_default_FSL" ]; then
    warn_msg
    warn_msg "[WARNING]: This version of CBIG repository's default FSL version is $CBIG_default_FSL."
    warn_msg "You appear to be using an unrecognized/invalid FSL version. \
Switch to FSL $CBIG_default_FSL if possible."
  elif [ "$USER_default_FSL" != "$CBIG_default_FSL" ]; then
    warn_msg
    warn_msg "[WARNING]: This version of CBIG repository's default FSL version is $CBIG_default_FSL."
    warn_msg "You appear to be using FSL $USER_default_FSL (which may or may not work with current repo). \
Switch to FSL version $CBIG_default_FSL if possible."
    warn_msg "Note: stable projects may not follow default setting, \
refer to the proper config file of the project you use."
    warn_msg "If you are unsure, please run our stable projects' examples \
to check if the results can be replicated."
  fi
fi

#############
# set up ANTS
#############

if [ -z $CBIG_ANTS_DIR ] || [ ! -d $CBIG_ANTS_DIR ]; then
  warn_msg
  warn_msg "[WARNING]: CBIG_ANTS_DIR is not set or points to a non-existent directory"
else
  export ANTSPATH=$CBIG_ANTS_DIR
  export PATH=$CBIG_ANTS_DIR:$PATH
fi

############################
# set up Workbench directory
############################

if [ -z $CBIG_WB_DIR ] || [ ! -d $CBIG_WB_DIR ]; then
  warn_msg
  warn_msg "[WARNING]: CBIG_WB_DIR is not set or points to a non-existent directory"
else
  export PATH=$CBIG_WB_DIR/bin_rh_linux64:$PATH

  # check whether user's default workbench version matches CBIG default workbench version
  CBIG_default_WB=`cat $CBIG_CODE_DIR/setup/CBIG_sample_config.sh | grep CBIG_WB_DIR \
    | grep -E -o [0-9]+.[0-9]+.[0-9]+`
  USER_default_WB=`wb_command | grep Version | grep -E -o [0-9]+.[0-9]+.[0-9]+`

  if [ -z "$USER_default_WB" ]; then
    warn_msg
    warn_msg "[WARNING]: This version of CBIG repository's default WORKBENCH version is $CBIG_default_WB."
    warn_msg "You appear to be using an unrecognized/invalid WORKBENCH version. \
Switch to WORKBENCH $CBIG_default_WB if possible."
  elif [ "$USER_default_WB" != "$CBIG_default_WB" ]; then
    warn_msg
    warn_msg "[WARNING]: This version of CBIG repository's default WORKBENCH version is $CBIG_default_WB."
    warn_msg "You appear to be using WORKBENCH $USER_default_WB (which may or may not work with current repo). \
Switch to WORKBENCH version $CBIG_default_WB if possible."
    warn_msg "Note: stable projects may not follow default setting, \
refer to the proper config file of the project you use."
    warn_msg "If you are unsure, please run our stable projects' examples \
to check if the results can be replicated."
  fi
fi

#############
# set up AFNI
#############

if [ -z $CBIG_AFNI_DIR ] || [ ! -d $CBIG_AFNI_DIR ]; then
  warn_msg
  warn_msg "[WARNING]: CBIG_AFNI_DIR is not set or points to a non-existent directory"
else
  export PATH=$CBIG_AFNI_DIR:$PATH

  # check whether user's default afni version matches CBIG default afni version
  CBIG_default_AFNI=`cat $CBIG_CODE_DIR/setup/CBIG_sample_config.sh | grep CBIG_AFNI_DIR \
    | grep -o afni.* | grep -o AFNI.* | sed 's/\/linux_openmp_64//g'`
  USER_default_AFNI=`afni --version | grep -o "AFNI_.*" | sed 's/)//g'`

  if [ -z "$USER_default_AFNI" ]; then
    warn_msg
    warn_msg "[WARNING]: This version of CBIG repository's default AFNI version is $CBIG_default_AFNI."
    warn_msg "You appear to be using an unrecognized/invalid AFNI version. \
Switch to AFNI $CBIG_default_AFNI if possible."
  elif [ "$USER_default_AFNI" != "$CBIG_default_AFNI" ]; then
    warn_msg
    warn_msg "[WARNING]: This version of CBIG repository's default AFNI version is $CBIG_default_AFNI."
    warn_msg "You appear to be using AFNI $USER_default_AFNI (which may or may not work with current repo). \
Switch to AFNI version $CBIG_default_AFNI if possible."
    warn_msg "Note: stable projects may not follow default setting, \
refer to the proper config file of the project you use."
    warn_msg "If you are unsure, please run our stable projects' examples \
to check if the results can be replicated."
  fi
fi

####################################
# create symlink to git hook scripts
####################################

if [ ! -f "$CBIG_CODE_DIR/.git/hooks/pre-commit" ];
then
  mkdir -p $CBIG_CODE_DIR/.git/hooks
  ln -s "$CBIG_CODE_DIR/hooks/pre-commit" "$CBIG_CODE_DIR/.git/hooks/pre-commit"
fi
if [ ! -f "$CBIG_CODE_DIR/.git/hooks/pre-push" ];
then
  mkdir -p $CBIG_CODE_DIR/.git/hooks
  ln -s "$CBIG_CODE_DIR/hooks/pre-push" "$CBIG_CODE_DIR/.git/hooks/pre-push"
fi

