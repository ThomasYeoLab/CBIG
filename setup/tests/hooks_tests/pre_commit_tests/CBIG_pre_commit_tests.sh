#!/bin/sh
curr_dir=`pwd`
for folder in */ 
do 
    echo "[TEST](start) $folder"
    cd $curr_dir
    rsync -az $folder $CBIG_CODE_DIR/stable_projects/
    git add $CBIG_CODE_DIR/stable_projects/$folder/*
    cd $CBIG_CODE_DIR
    sh $CBIG_CODE_DIR/hooks/pre-commit
    git reset
    rm -r $CBIG_CODE_DIR/stable_projects/$folder
    echo "You should get following output:"
    case $folder in
        A_check_CBIG_prefix_scripts/)
            echo "   [FAILED] There are functions without CBIG_ prefix. Abort committing."
            ;;
        B_check_MIT_license/)
            echo "   [FAILED] There are functions without or not following our MIT license. Abort committing."
            ;;
        C_check_addpath_rmpath/)
            echo "   [FAILED] There are functions without 'rmpath' at the end. Abort committing."
            ;;
        D_check_CBIG_prefix_matlab_class/)
            echo "   [FAILED] There are Matlab classes without CBIG_ prefix. Abort committing."
            ;;
        E_check_character_limit/)
            echo "   [FAILED] There are scripts which does not follow 120 characters limit of each line. Abort committing."
            ;;
        F_check_pep8_format/)
            echo "   [FAILED] There are python scripts which does not follow yapf pep8 format. Abort committing."
    esac
    echo "[TEST](end) $folder"
    echo ""
    read -p "Press enter to continue"
done 
