## Background

This folder contains the MRtrix scripts to run using CSCHPC's cluster. For full details on MRtrix3, please visit their website
[https://mrtrix.readthedocs.io/en/latest/](https://mrtrix.readthedocs.io/en/latest/).

### Usage

The MRtrix pipeline compromises of 2 main steps. 

- Preparation prior to tratography

  Please install MRtrix in your python environment prior to running these scripts. You can either install the python environment with `MRtrix_env.txt` 
  in this folder or the CBIG python 3 environment found under `$CBIG_CODE_DIR/setup/python_env_setup/CBIG_python3_env.yml` (recommended). 

  Prepare the data by creating a folder titled `diffusion` containing 
  each subject's diffusion image, bval and bvec file. They should all have the same root name. Run Freesurfer's `recon_all` on each subject's T1 image 
  and place them in a `recon_all` folder (creating symbolic links from the original recon_all folder will save you the trouble of moving them). 
  Additionally, if you have a path of b=0 mask files, you can specify it in order to skip the b=0 mask calculation step. Each subject should have it's 
  own folder in the mask directory and be named as ${subject_id}_bet_b0_mask.nii.gz. The file structure for the input directory should be as follows: 

```
Input directory (input_dir)
|-- diffusion
    |-- sub-1
    |-- sub-2
       |-- sub-2.nii.gz
       |-- sub-2.bval
       |-- sub-2.bvec
|-- recon_all
    |-- fsaverage
    |-- sub-1
    |-- sub-2 
       |-- label
       |-- mri
       |-- mri.2mm
       |-- scripts
       |-- surf

Mask directory (mask_output_dir)
|-- sub-1
|-- sub-2
|-- sub-3
    |-- sub-1_bet_b0_mask.nii.gz

```

- Step 2: Running the script

  Prepare a list of subjects to run. The script can then be run with the following command. Remember to specify the same python environment as 
  your installation:
```
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG2022_DiffProc/ \
		MRtrix/CBIG_DiffProc_batch_tractography.sh \
		--subj_list /path/to/txtfile --dwi_dir /path/to/dwi_images \
		--output_dir /path/to/output --py_env name_of_AMICO_environment \
		--mask_output_dir /path/to/stored/b0_brainmask
```

   In the first step, `CBIG_DiffProc_preprocess_tractography` converts each parcellation to the native volume space in nifti format, and then converts 
   it to the mif format. By default, it retrieves the Schaefer 100 and 400 parcellations, as well as the Desikan-Killiany parcellation. Tracts are 
   generated for the diffusion data. You will need the annot files for the Schaefer parcellations, and the Freesurfer recon-all have been run already to
   generate the Desikan-Killiany parcellation.

   By default, the iFOD2 probablistic algorithm is used. If it is preferred, the user can select to do deterministic 
   tractography by specifying `algo` to be FACT. During tract generation, we use anatomy constrained tractography using the freesurfer segmentation of the 
   white and gray matter. A detailed list of settings for tract generation can be found in `CBIG_DiffProc_tractography_2_gen_tractogram.sh` under the 
   `tckgen` commands.

   In the second step, `CBIG_DiffProc_tractography_3_gen_connectome` samples each tract and generates the the structural connectivity stream count.
   Additionally, a DTI model is fit to the model and additional metrics can be sampled along each tract. The script needs to be modified if other 
   metrics (e.g. from NODDI) are to be included. After this, `CBIG_DiffProc_mrtrix_fill_na.m` replaces any missing parcels with NaN.

   Lastly, `CBIG_DiffProc_tractography_4_del_tractograms` can be used to clean up all files except the connectome if the files are too large.

----

## Bugs and Questions

Please contact Leon Ooi at leonooiqr@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.