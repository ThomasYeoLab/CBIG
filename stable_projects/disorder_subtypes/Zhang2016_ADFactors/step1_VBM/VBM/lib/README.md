# FSL-VBM -- General

**Note** This is the general procedure, which is easier to be applied to a new dataset. To replicate our PNAS results or view a simplified version for concept/procedure understanding, please refer to folder `../replicatePNAS/` instead.

Given brain images, this folder runs FSL-VBM, as described at http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLVBM/UserGuide.

---
## Master Scripts `CBIG_runVBM.sh` and `CBIG_runVBM_givenTmp.sh`

`CBIG_runVBM.sh` is the master script that runs each step of FSL-VBM, whereas `CBIG_runVBM_givenTmp.sh` skips the VBM steps that construct the study-specific template (i.e., runs Steps 1, 6, 7 and 8). 

If you have a cluster, both scripts will parallelize the jobs whenever it is possible. 

`../replicatePNAS/CBIG_runVBM_wrapper.sh` demonstrates how to call `CBIG_runVBM.sh` and `CBIG_runVBM_givenTmp.sh`. You may also do `./CBIG_runVBM.sh` (or `./CBIG_runVBM_givenTmp.sh`) for detailed usage.

----
## Step 1: GM Segmentation

`CBIG_step1_segGM.sh` is the first script that `CBIG_runVBM.sh` calls. This script either runs `CBIG_step1_segGM_job.sh` serially or submits it to the job queue if you have a cluster.

At the end of this step, gray matter (GM) will have been segmented out for each brain image. The output files (GM images) are identified as `${ID}_pve_1.nii.gz`.

In the case where GM segmentation is the only step that you want to run, type `./CBIG_step1_segGM.sh` to see detailed usage of this script. If you are unsatisfied with the FAST (the FSL function for segmentation) parameters, edit `./CBIG_step1_segGM_job.sh`.

----
## Step 2: Affine Registration to Standard Template

`CBIG_step2_affineRegToStdTmp.sh` is the second script that `CBIG_runVBM.sh` calls. This script either runs `CBIG_step2_affineRegToStdTmp_job.sh` serially or submits it to the job queue if you have a cluster.

At the end of this step, GM images selected for building the study-specific template will have been registered to the standard FSL GM template (or any other standard template that you provide). The output files (GM images registered to a standard template) are identified as `${ID}_GMToStdTmp.nii.gz`.

In the case where affine registration is the only step that you want to run, type `./CBIG_step2_affineRegToStdTmp.sh` to see detailed usage of this script.

----
## Step 3: Creating Study-Specific Affine Template

`CBIG_step3_createAffineTmp.sh` is the third script that `CBIG_runVBM.sh` calls.

At the end of this step, a study-specific, left-right-symmetric affine template will be created with the registered GM images from the previous step (`*_GMToStdTmp.nii.gz`). The output file (affine template) is identified as `affineTmp.nii.gz`.

In the case where creating an affine template is the only step that you want to run, type `./CBIG_step3_createAffineTmp.sh` to see detailed usage of this script. 

----
## Step 4: Nonlinear Registration to the Affine Template

`CBIG_step4_nonlinRegToAffineTmp.sh` is the fourth script that `CBIG_runVBM.sh` calls. This script either runs `CBIG_step4_nonlinRegToAffineTmp_job.sh` serially or submits it to the job queue if you have a cluster.

At the end of this step, the same set of GM images (those selected for building the study-specific template) will have been nonlinearly registered to the affine template from the previous step (i.e., `affineTmp.nii.gz`). The output files (GM images registered to the affine template) are identified as `${ID}_GMToAffineTmp.nii.gz`.

In the case where nonlinear registration is the only step that you want to run, type `./CBIG_step4_nonlinRegToAffineTmp.sh` to see detailed usage of this script.

----
## Step 5: Creating Study-Specific Nonlinear Template

`CBIG_step5_createNonlinTmp.sh` is the fifth script that `CBIG_runVBM.sh` calls.

At the end of this step, a study-specific, left-right-symmetric nonlinear template will be created with the registered GM images from the previous step (`*_GMToAffineTmp.nii.gz`). The output file (nonlinear template) is identified as `nonlinTmp.nii.gz`.

In the case where creating a nonlinear template is the only step that you want to run, type `./CBIG_step5_createNonlinTmp.sh` to see detailed usage of this script.

----
## Step 6: Nonlinear Registration to the Nonlinear Template (for *all* GM images)

`CBIG_step6_nonlinRegToNonlinTmp.sh` is the sixth script that `CBIG_runVBM.sh` calls. This script either runs `CBIG_step6_nonlinRegToNonlinTmp_job.sh` serially or submits it to the job queue if you have a cluster.

At the end of this step, *all* GM images (rather than only the GM images selected for building the template) will have been nonlinearly registered to the nonlinear template from the previous step (i.e., `nonlinTmp.nii.gz`). The output files (GM images registered to the nonlinear template) are identified as `${ID}_GMToNonlinTmp_mod.nii.gz`, where `_mod` indicates modulation due to the nonlinear registration. For more details, see the FSL-VBM user guide.

In the case where nonlinear registration is the only step that you want to run, type `./CBIG_step6_nonlinRegToNonlinTmp.sh` to see detailed usage of this script.

----
## Step 7: Concatenating and Generating a Binary GM Mask

`CBIG_step7_concatAndGenMask.sh` is the seventh script that `CBIG_runVBM.sh` calls.

At the end of this step, a 4D volume (identified as `GMToNonlinTmp_mod_4d.nii.gz`) will be produced by concatenating the registered (and modulated) GM images from the previous step (`*_GMToNonlinTmp_mod.nii.gz`). In addition, a binary GM mask (identified as `GMToNonlinTmp_mod_mean_binThr0.05.nii.gz`) will be created by averaging that 4D volume and then binarizing at 0.05.

In the case where this is the only step that you want to run, type `./CBIG_step7_concatAndGenMask.sh` to see detailed usage of this script. If you want to change the binarization threshold from 0.05 to some other values, edit `CBIG_step7_concatAndGenMask.sh`.

----
## Step 8: Smoothing

`CBIG_step8_smooth.sh` is the final script that `CBIG_runVBM.sh` calls.

At the end of this step, each registered GM image (`${ID}_GMToNonlinTmp_mod.nii.gz`) will be smoothed by a Gaussian kernel whose sigma can be specified. The output file (smoothed images) is identified as `GMToNonlinTmp_mod_4d_s${sigma}.nii.gz`.

In the case where smoothing is the only step that you want to run, type `./CBIG_step8_smooth.sh` to see detailed usage of this script.

----
## Helper Function `CBIG_waitUntilFinished.sh`

This helper function is useful when you use a cluster to parallelize your jobs. For example, the completion of Step 2 is a prerequisite for Step 3, and jobs in Step 2 are run in parallel as jobs in the cluster. In this case, this helper function holds off Step 3 until all jobs in Step 2 are finished.
