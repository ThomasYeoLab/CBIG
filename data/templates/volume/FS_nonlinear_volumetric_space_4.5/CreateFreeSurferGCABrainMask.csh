#! /bin/csh -f

# setting up parameters
#set SUBJECTS_DIR = $NEXUS_GSP
#setenv SUBJECTS_DIR $NEXUS_GSP


# Create atlas labels and mean
mri_convert $FREESURFER_HOME/average/RB_all_2008-03-26.gca#1 gca_labels2mm.mgz
mri_convert $FREESURFER_HOME/average/RB_all_2008-03-26.gca#0 gca_mean2mm.mgz

# Create Thalamus Mask
mri_binarize --i gca_labels2mm.mgz --match 9 --match 10 --match 49 --match 48 --o ThalamusMask.GCA.t0.5.nii.gz

# Create Cerebellum Mask
mri_binarize --i gca_labels2mm.mgz --match 8 --match 47 --o CerebellumMask.GCA.t0.5.nii.gz
mri_binarize --i gca_labels2mm.mgz --match 7 --match 8 --o lh.CerebellumMask.GCA.t0.5.nii.gz
mri_binarize --i gca_labels2mm.mgz --match 46 --match 47 --o rh.CerebellumMask.GCA.t0.5.nii.gz

# Create Basal Ganglia Mask
mri_binarize --i gca_labels2mm.mgz --match 11 --match 12 --match 13 --match 26 --match 27 --match 50 --match 51 --match 52  --match 58 --match 59 --o BasalGangliaMask.GCA.t0.5.nii.gz

# Create Hippocampus Mask
mri_binarize --i gca_labels2mm.mgz --match 17 --match 53 --o HippocampusMask.GCA.t0.5.nii.gz

# Create Amygdala Mask
mri_binarize --i gca_labels2mm.mgz --match 18 --match 54 --o AmygdalaMask.GCA.t0.5.nii.gz

# Create VentralDC Mask
mri_binarize --i gca_labels2mm.mgz --match 28 --match 60 --o VentralDCMask.GCA.t0.5.nii.gz

# Create BrainStem Mask
mri_binarize --i gca_labels2mm.mgz --match 16 --match 15 --o BrainStemMask.GCA.t0.5.nii.gz

# Create Striatum Mask
mri_binarize --i gca_labels2mm.mgz --match 11 --match 12 --match 26 --match 50 --match 51 --match 58 --o StriatumMask.GCA.t0.5.nii.gz
mri_binarize --i gca_labels2mm.mgz --match 11 --match 12 --match 26 --o lh.StriatumMask.GCA.t0.5.nii.gz
mri_binarize --i gca_labels2mm.mgz --match 50 --match 51 --match 58 --o rh.StriatumMask.GCA.t0.5.nii.gz

# Create Caudate Mask
mri_binarize --i gca_labels2mm.mgz --match 11 --o lh.caudate.GCA.t0.5.nii.gz
mri_binarize --i gca_labels2mm.mgz --match 50 --o rh.caudate.GCA.t0.5.nii.gz

# Create Putamen Mask
mri_binarize --i gca_labels2mm.mgz --match 12 --o lh.putamen.GCA.t0.5.nii.gz
mri_binarize --i gca_labels2mm.mgz --match 51 --o rh.putamen.GCA.t0.5.nii.gz


# Create Loose GCA head mask
set cmd = "matlab -nojvm -nodesktop -nosplash -r 'CBIG_CreateLooseGCAHeadMask gca_labels2mm.nii.gz LooseHeadMask.GCA.t0.5.nii.gz 2'"
echo $cmd
eval $cmd

# Create Loose GCA head mask
set cmd = "matlab -nojvm -nodesktop -nosplash -r 'CBIG_CreateLooseGCAHeadMask gca_labels2mm.nii.gz ReallyLooseHeadMask.GCA.t0.5.nii.gz 4'"
echo $cmd
eval $cmd

# Mask dilation
set cmd = "matlab -nojvm -nodesktop -nosplash -r 'CreateLooseSubcorticalMasks'"
echo $cmd
eval $cmd

# transform masks 
foreach mask (LooseCerebellumMask.GCA.t0.5 LooseCerebellum.dist1.Mask.GCA.t0.5)
    TransformMNI2mm2FS2mm.csh $mask.nii.gz $mask.MNI152_2mm.nii.gz nearest 0
end


# Create Brain Mask
foreach threshold (0.5)

    set input_text_file = $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/input.txt
    rm $input_text_file;

    echo "$CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/gca_labels2mm.mgz" >> $input_text_file
	    
    set output_file = $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/BrainMask.GCA.t$threshold.nii.gz   
    set cmd = "matlab -nojvm -nodesktop -nosplash -r 'CreateFreeSurferBrainMask $output_file $threshold $input_text_file'"
    echo $cmd
    eval $cmd
    rm $input_text_file
end

# Create Subcort Mask
foreach threshold (0.5)

    set input_text_file = $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/input.txt
    rm $input_text_file;

    echo "$CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/gca_labels2mm.mgz" >> $input_text_file
        
    set output_file = $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/SubcortMask.GCA.t$threshold.nii.gz   
    set cmd = "matlab -nojvm -nodesktop -nosplash -r 'CreateFreeSurferSubcortMask $output_file $threshold $input_text_file'"
    echo $cmd
    eval $cmd
    rm $input_text_file
end

foreach threshold (0.5)

    set input_text_file = $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/input.txt
    rm $input_text_file;

    echo "$CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/gca_labels2mm.mgz" >> $input_text_file
        
    set output_file = $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/SubcortMask.GCA.t$threshold.nii.gz   
    set cmd = "matlab -nojvm -nodesktop -nosplash -r 'CreateFreeSurferSubcortMask $output_file $threshold $input_text_file'"
    echo $cmd
    eval $cmd
    rm $input_text_file
end

# Create Brain + Subcort Mask: Subcort = 1, Rest of Brain = 2
fslmaths $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/BrainMask.GCA.t$threshold.nii.gz -add $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/SubcortMask.GCA.t$threshold.nii.gz BrainSubcortMask.GCA.t$threshold.nii.gz


# Subcort Mask + cerebellum white matter
foreach threshold (0.5)

    set input_text_file = $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/input.txt
    rm $input_text_file;

    echo "$CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/gca_labels2mm.mgz" >> $input_text_file
        
    set output_file = $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/SubcortCerebellumWhiteMask.GCA.t$threshold.nii.gz   
    set cmd = "matlab -nojvm -nodesktop -nosplash -r 'CreateFreeSurferSubcortAndCerebellumWhiteMask $output_file $threshold $input_text_file'"
    echo $cmd
    eval $cmd
    rm $input_text_file
end

# Subcort Mask + cerebellum white matter = 1 ; Rest of brain = 2 
fslmaths $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/BrainMask.GCA.t$threshold.nii.gz -add $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/SubcortCerebellumWhiteMask.GCA.t$threshold.nii.gz BrainSubcortCerebellumWhiteMask.GCA.t$threshold.nii.gz


# Subcort Mask + but now with Loose Striatum Mask added
fslmaths $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/SubcortMask.GCA.t$threshold.nii.gz -add $CODE_DIR/templates/volume/FS_nonlinear_volumetric_space_4.5/LooseStriatum.dist1.Mask.GCA.t0.5.nii.gz SubcortLooseStriatum.dist1.Mask.GCA.t0.5.nii.gz

mri_binarize --i SubcortLooseStriatum.dist1.Mask.GCA.t0.5.nii.gz --match 1 --match 2 --o SubcortLooseStriatum.dist1.Mask.GCA.t0.5.nii.gz
