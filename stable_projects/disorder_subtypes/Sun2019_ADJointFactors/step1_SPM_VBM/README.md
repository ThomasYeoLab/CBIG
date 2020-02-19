## SPM VBM
SPM stands for Statistical Parametric Mapping which is a popular neuroscitific software. VBM stands for voxel-based morphometry which is a neuroimaging analysis technique that allows investigation of focal differences in brain anatomy. If you are not familar with SPM VBM, please read the `doc/SPM_VBM_pipeline.pdf` to understand SPM VBM workflow. 

----
## Code
All the code of SPM VBM are in `code` folder. The main functions of SPM VBM are `code/CBIG_MMLDA_runVBM.sh`, which will create a new study specific template and `code/CBIG_MMLDA_runVBM_givenTemp.sh`, which will use existing template. The following list all steps used in `code/CBIG_MMLDA_runVBM.sh`. For `code/CBIG_MMLDA_runVBM_givenTemp.sh`, we don't need to run step2-5. 

1. Convert raw image files to .nii files for processing
  * Inputs
    *  List of subject IDs (.txt)
    *  List of paths to subject raw files (.txt)
    *  Output directory
  * Output
    *  SubjectID.nii files in output directory

1a. Reorientation of image files so that the position of origin is in AC (anterior commissure) [Compulsory] 
  * Input
    *  Subject files (xxx.nii)
    *  Sample reoriented file (xxx.mat)
  * Output
    *  Reoriented subject files (SubjectID_reorient.nii)
    *  Reorientation matrix (SubjectID_reorient.mat)
  * Orientation QC
    *  Use spm GUI and display each MRI scan as mentioned in the following Youtube video.
  * Manual Steps
    1. Follow the [youtube video](https://www.youtube.com/watch?v=AwNJAUKLhqY) to reorient one sample image and you will get xxx.mat which includes the [translation matrix](https://en.wikipedia.org/wiki/Translation_%28geometry%29).
    2. By using step1a, you will first do a coarse transformation by apply the translation matrix to the image and then do a fine transformation which is [rigid transformation](https://en.wikipedia.org/wiki/Rigid_transformation) from the image to template. These two steps will make sure the origin is near the AC.
    When you run step1a at first time, you can try `Translation_Indentity.mat` which will not do any translation but just do rigid tranformation.
    3. After applying step1a once, a batch of images may succeed and it will output `SubjectID_reorient.mat`. For the remaining failed subjects, you have to repeat 1~2 until all subjects have the `SubjectID_reorient.mat`.
    4. Now you need to change your subject list by adding `_reorient` after each Subject ID and then run remaining steps.
    PS: If you encounter the "Insufficient image overlap" error when doing step2 segmentation, it means that the orientation of the T1 image is wrong. Then, you need to go back and check why it is wrong and reorient it.  

2. First segmentation using default templates
  * Input
    *  Subject files (.nii)
    *  Tissue probability map (tpm.nii) from SPM
    *  DARTEL tissue probability map (Template_1_IXI555_MNI152.nii) from SPM toolbox
  * Output
    *  Grey matter, white matter, and cerebrospinal fluid segments in native space (p1xxxxx.nii, p2xxxxx.nii, p3xxxxx.nii)
    *  Grey matter, white matter, and cerebrospinal fluid segments after rigid registration to template (rp1xxxxx.nii, rp2xxxxx.nii, rp3xxxxx.nii)
    *  Grey matter, white matter, and cerebrospinal fluid segments after affine registration to template (rp1xxxxx_affine.nii, rp2xxxxx_affine.nii, rp3xxxxx_affine.nii)

2a. Segmentation QC
    *  The script 2a will call fsl slicesdir function to overlay segmentation GM to original T1 image, you can check the index.html in slicesdir folder.
    *  There are several parts you can check: 1) the red outline shouldn't include neck and eyes, 2) the red outline shouldn't miss any part of the brain, such as cerebellum, frontal lobe. 

3. Creating customized DARTEL template for resegmentation
  * Input
    *  Grey and white matter segments after affine registration (rp1xxxxx_affine.nii, rp2xxxxx_affine.nii)
  * Output
    *  (u_rp1xxxxxaffine_Template.nii)
    *  Iterative steps of DARTEL normalization (Template_0-6.nii)

4. Population to ICBM registration
  * Input
    *  Final DARTEL Template (Template_6.nii)
  * Output
    *  DARTEL Template after ICBM registration (y_Template_6_2mni.nii)
 
5. Deformation
  * Input
    *  DARTEL Template after ICBM registration (y_Template_6_2mni.nii)
    *  Iterative templates of DARTEL normalization (Template_0-6.nii)
  * Output
    *  Templates after deformation (wTemplate_0-6.nii)

6. Resegmentation using customized DARTEL template
  * Input
    *  Subject files (xxxxx.nii)
    *  Tissue probability map (tpm.nii) from SPM toolbox
    *  Customized DARTEL tissue probability map (wTemplate_1.nii)
  * Output
    *  Grey matter after resegmentation (mwp1xxxxx.nii)
    *  Intermediate output (wmxxxxx.nii)

7. Smoothing
  * Input
    *  Grey matter after resegmentation (mwp1xxxxx.nii)
  * Output
    *  Smoothed results (s10_mwp1xxxxx.nii)

8. Merge and generate mask
  * Input
    *  Original (mwp1xxxxx.nii) and smoothed outputs (s10_mwp1xxxxx.nii)
  * Output
    *  Merged original (gm_merg.nii.gz) and smoothed outputs (gm_merg_s10.nii.gz)
    *  Mask (mask_fsl_0_1.nii.gz)

9. Calculate Grey Matter and Intracranial Volume
  * Input
    *  Grey matter, white matter, and cerebrospinal fluid segments in native space (p1xxxxx.nii, p2xxxxx.nii, p3xxxxx.nii)
  * Output
    *  Matlab file with GM and ICV values (id_gmVol_icv.csv)

10. Downsample from 1.5mm to 2mm
  * Input
    *  Grey matter segments in MNI space (mwp1xxx.nii)
  * Output
    *  Grey matter segments in MNI 2mm space (mwp1xxx_MNI2mm.nii)

11. Smoothing the grey matter segemnts in 2mm space
  * Input
    *  Grey matter segments in MNI 2mm space (mwp1xxx_MNI2mm.nii) 
  * Output
    *  Smoothed results (s10_mwp1xxx_MNI2mm.nii)

12. Merge and generate mask
  * Input
    *  Original (mwp1xxx_MNI2mm.nii) and smoothed outputs (s10_mwp1xxx_MNI2mm.nii)
  * Output
    *  Merged original (gm_merg_MNI2mm.nii.gz) and smoothed outputs(gm_merg_MNI2mm_s10.nii.gz)
    *  Mask (mask_fsl_0_1_MNI2mm_maskbrain.nii.gz)
----
## Examples
The examples are here: `code/CBIG_MMLDA_runVBM_example.sh` `code/CBIG_MMLDA_runVBM_givenTemp_example.sh`. You can also find the results of examples in `${CBIG_TESTDATA_DIR}/stable_projects/disorder_subtypes/Sun2018_ADJointFactors/step1_SPM_VBM/results`.

