#!/bin/sh
curr_dir=`pwd`
for folder in */ 
do 
    echo "[TEST](start) $folder"
    cd $curr_dir
    rsync -az $folder $CBIG_CODE_DIR/stable_projects/
    git add $CBIG_CODE_DIR/stable_projects/$folder/*
    git commit -m "test"
    cd $CBIG_CODE_DIR
    yes | $CBIG_CODE_DIR/hooks/pre-push
    git reset HEAD~1
    rm -r $CBIG_CODE_DIR/stable_projects/$folder
    echo "[TEST](end) $folder"
    echo ""
    read -p "Press enter to continue"
done 
