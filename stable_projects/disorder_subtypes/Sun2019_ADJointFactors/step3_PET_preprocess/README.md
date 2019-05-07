## PET preprocessing
PET stands for positron-emission tomography is a nuclear medicine functional imaging technique that is used to observe metabolic processes in the body. In our paper, we use the tau PET imaging with 18F-AV1451 as the tracer to do analyses. 

----

## Code 
The main function of PET preprocessing is `CBIG_MMLDA_runPETpreprocess.sh`. The following list all steps within the main function.

1. Convert raw T1 image and PET image to nifti file and rename it
    - Input:
        - List of subjects IDs (.txt)
        - List of paths to T1 image of corresponding subjects (.txt)
        - List of paths to PET image of corresponding subjects (.txt)
        - Output directory
    - Output:
        - T1 image: `${out_dir}/${id}/T1.nii.gz`
        - PET image:`${out_dir}/${id}/PET.nii.gz`

2. Freesurfer recon-all of T1 image
    - Input:
        - List of subjects IDs (.txt)
        - Output directory
    - Output:
        - recon-all folder: `${out_dir}/${id}_FS`

3. Freesurfer register PET to T1
    - Input:
        - List of subjects IDs (.txt)
        - Output directory
    - Output:
        - PET to T1 registration file: `${out_dir}/${id}/PET2T1.reg.lta`
        - PET in T1 space: `${out_dir}/${id}/PET_in_T1.nii.gz`

4. Normalize PET image wrt cerebellum gm
    - Input:
        - List of subjects IDs (.txt)
        - Script directory
        - Output directory
    - Output:
        - Normalize PET with respect to cerebellum grey matter: `${out_dir}/${id}/PET_in_T1_NormCere.nii.gz`

5. Transform PET image in T1 space to MNI space
    - Input:
        - List of subjects IDs (.txt)
        - Directory of deformation field
        - Script directory
        - SPM directory
        - Output directory
    - Output:
        - Reorient the image based on T1: `${out_dir}/${id}/PET_in_T1_NormCere_reorientT1.nii.gz`
        - Transfrom the image to MNI space: `${out_dir}/${id}/wPET_in_T1_NormCere_reorientT1.nii`
        - Downsample to MNI2mm: `${out_dir}/${id}/PET_in_T1_NormCere_SPMVBMdeform_MNI2mm.nii.gz`

6. Merge PET image of subjects
    - Input:
        - List of subjects IDs (.txt)
        - Script directory
        - Output directory
    - Output:
        - Merged file: `${out_dir}/merg_SPMVBMdeform.nii.gz`

----

## Examples
The examples script is here: `CBIG_MMLDA_runPETpreprocess_example.sh`. For users outside CBIG lab, you may not be able to run it because we use ADNI data here. However, you can use it as reference. For users inside CBIG lab, you can just run it and see whether it [PASSED] or [FAILED] and it will take about 12hours to run due to the recon-all.
