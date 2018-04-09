Reference
=========
Yeo BT, Krienen FM, Sepulcre J, Sabuncu MR, Lashkari D, Hollinshead M, Roffman JL, Smoller JW, Zollei L., Polimeni JR, Fischl B, Liu H, Buckner RL. [**The organization of the human cerebral cortex estimated by intrinsic functional connectivity**](http://www.ncbi.nlm.nih.gov/pubmed/21653723). Journal of Neurophysiology 106(3):1125-1165, 2011.

----

Background
==========
Resting state fMRI data from 1000 subjects were registered using surface-based alignment. A clustering approach was employed to identify and replicate 7 and 17 networks of functionally coupled regions across the cerebral cortex. The results revealed local networks confined to sensory and motor cortices as well as distributed networks of association regions that form interdigitated circuits. Within the sensory and motor cortices, functional connectivity followed topographic representations across adjacent areas. In association cortex, the connectivity patterns often showed abrupt transitions between network boundaries, forming largely parallel circuits.

----

Information about Downloads
===========================
There are three folders in `Yeo_JNeurophysiol11_SplitLabels`, corresponding to the "fsaverage5" surface space, "MNI152" volumetric space and scripts directory whose scripts were used to generate the MNI152 parcellations from the surface parcellations.

----

Information about fsaverage5 folder
===================================
1) The structure of fsaverage5 follows that of a preprocessed freesurfer subject. In particular, `fsaverage5/label/` contain all the parcellation and confidence files. For example, `fsaverage5/label/rh.Yeo2011_7Networks_N1000.annot` is the 7-network parcellation for 1000 subjects on the right hemisphere and `fsaverage5/label/lh.Yeo2011_17NetworksConfidence_N1000.mgz` is the confidence map for the 17-network parcellation for 1000 subjects on the left hemisphere. 

2) The `*.split_components.annot` are obtained by

   a) performing a connected component analysis of the original parcellation
   
   b) for components containing <= 7 vertices, vertices are reassigned to most closely correlated networks
   
   c) certain large components in the 17 networks are further split. For example, since we have highly accurate probability maps of histological V1, we use it to split the visual networks into striate and extra-striate components 
   
   d) vertices at boundary (between networks) are peeled back

3) The labels corresponding to the `*.split_components.annot` are found in `fsaverage5/label/split_labels*/`

----

Example Usage of fsaverage5
===========================
1) Make sure freesurfer has been set up and configured as usual (http://surfer.nmr.mgh.harvard.edu/fswiki/SetupConfiguration).

2) In terminal,

   a) "cd" to unzipped folder containing this README file, fsaverage5
   
   b) tksurfer fsaverage5 lh inflated -annotation fsaverage5/label/lh.Yeo2011_7Networks_N1000.split_components.annot

3) The split labels are not symmetric, so that for example, lh.7Networks_LH_Cont_PFCl.label is not necessary the left hemisphere homologue of rh.7Networks_LH_Cont_PFCl.label

----

Information about MNI152 folder
===============================
1) **Yeo2011_??Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii.gz** is the projection of the surface-based split annotation into `norm.nii.gz` under `$CBIG_CODE_DIR/CBIG_private/data/templates/volume/FSL_MNI152_FS4.5.0/mri`, which is the FSL MNI152 1mm template interpolated and intensity normalized into a 256 x 256 x 256 1mm-isotropic volume (obtained by putting the FSL MNI152 1mm template through recon-all using FreeSurfer 4.5.0).

2) **Yeo2011_??Networks_N1000.split_components.FSL_MNI152_2mm.nii.gz** is the downsampling of Yeo2011_??Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii.gz into FSL MNI152 2mm space using template `MNI152_T1_2mm_brain.nii.gz` under `$CBIG_CODE_DIR/data/templates/volume/FSL5.0.8_MNI_templates`. Note that Yeo2011_??Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii.gz and Yeo2011_??Networks_N1000.split_components.FSL_MNI152_2mm.nii.gz are in the same space, but have different resolution.

3) **Yeo2011_??Networks_N1000.split_components.FSL_MNI152_1mm.nii.gz** is the resampling of Yeo2011_??Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii.gz into FSL MNI152 1mm space using template `MNI152_T1_1mm_brain.nii.gz` under `$CBIG_CODE_DIR/data/templates/volume/FSL5.0.8_MNI_templates`.

4) **??Networks_ColorLUT_freeview.txt** are colortables that can be used with FreeSurfer viewer freeview.

5) **??Networks_ColorLUT_fslview.lut** are colortables that can be used with FSL fslview.

----

Example Usage of MNI152
=======================
1) Except for the colortables, all the volumes are nifti volumes which can be read using any software like freeview (FreeSurfer), fslview (FSL), etc.

2) To overlay the 7-network volume on the templates in freeview (assuming FreeSurfer is already installed) with the appropriate colors, launch freeview with the following arguments, assuming the working directory is in the same directory as this README:
    ```
    freeview -v $CBIG_CODE_DIR/data/templates/volume/FSL5.0.8_MNI_templates/MNI152_T1_1mm_brain.nii.gz MNI152/Yeo2011_7Networks_N1000.split_components.FSL_MNI152_1mm.nii.gz:colormap=lut:lut=MNI152/7Networks_ColorLUT_freeview.txt
    ```
    or

    ```
    freeview -v $CBIG_CODE_DIR/data/templates/volume/FSL5.0.8_MNI_templates/MNI152_T1_2mm_brain.nii.gz MNI152/Yeo2011_7Networks_N1000.split_components.FSL_MNI152_2mm.nii.gz:colormap=lut:lut=MNI152/7Networks_ColorLUT_freeview.txt
    ```
    or
    
    ```
    freeview -v $CBIG_CODE_DIR/CBIG_private/data/templates/volume/FSL_MNI152_FS4.5.0/mri/norm.nii.gz MNI152/Yeo2011_7Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii.gz:colormap=lut:lut=MNI152/7Networks_ColorLUT_freeview.txt
    ```

3) To overlay the 7-network volume on the FSL MNI152 2mm template in fslview (assuming FSL is already installed) with the appropriate colors, assuming the working directory is in the same directory as this README:

   a) type `fslview $CBIG_CODE_DIR/data/templates/volume/FSL5.0.8_MNI_templates/MNI152_T1_2mm_brain.nii.gz MNI152/Yeo2011_7Networks_N1000.split_components.FSL_MNI152_2mm.nii.gz`
   
   b) In fslview, load the `7Networks_ColorLUT_fslview.lut` lookup table
   
   c) **IMPORTANT**: In the fslview window, set the minimum value of the volume to be 0 and maximum value to be 52. For 17 networks, the minimum value should be set to 0 and the maximum value to be set to 115

4) The names of the labels can be found in `*Networks_ColorLUT_freeview.txt`. Note that the labels are not symmetric, so that for example, 7Networks_LH_Cont_PFCl is not necessary the left hemisphere homologue of 7Networks_RH_Cont_PFCl.


5) To check the meaning of the name of the labels in `*Networks_ColorLUT_freeview.txt`, we provide two csv files under Yeo_JNeurophysiol11_SplitLabels: **Yeo2011_7networks_N1000.split_components.glossary.csv** and **Yeo2011_17networks_N1000.split_components.glossary.csv**.

----

Other Downloads
===============

1) Seed-based fcMRI Movies can be downloaded here: http://www.youtube.com/YeoKrienen

2) Cerebellar parcellation based on functional connectivity to cortex can be found here: http://surfer.nmr.mgh.harvard.edu/fswiki/CerebellumParcellation_Buckner2011

3) Seed-based fcMRI movies of the cerebellum with the cerebral cortex can be found here: http://www.youtube.com/bucknerkrienen

