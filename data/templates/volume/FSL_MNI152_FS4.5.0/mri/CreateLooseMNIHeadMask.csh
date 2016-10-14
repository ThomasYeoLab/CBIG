#! /bin/csh -f

# Create Loose GCA head mask
set cmd = "matlab -nojvm -nodesktop -nosplash -r 'CBIG_CreateLooseGCAHeadMask norm2mm.nii.gz LooseMNIHeadMask.2mm.nii.gz 2'"
echo $cmd
eval $cmd

# Create Loose GCA head mask
set cmd = "matlab -nojvm -nodesktop -nosplash -r 'CBIG_CreateLooseGCAHeadMask norm2mm.nii.gz ReallyLooseMNIHeadMask.2mm.nii.gz 4'"
echo $cmd
eval $cmd
