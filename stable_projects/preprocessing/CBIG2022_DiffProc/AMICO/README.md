## Background

This folder contains the scripts to run AMICO using CSCHPC's cluster. For full details on AMICO, please refer to Daducci et al., 2015 
[Accelerated Microstructure Imaging via Convex Optimization (AMICO) from diffusion MRI data](https://doi.org/10.1016/j.neuroimage.2014.10.026).

In brief, these scripts fit the NODDI model to multishell diffusion images. This estimates the intracellular and isotropic volume fractions, as well as the orientation dispersion at a voxel level. This pipeline creates 3 images of the aforementioned indices, which can then be used in other processes such as TBSS.

### Usage

Please install the AMICO package in your python environment prior to running these scripts. You can either install the python environment with 
`AMICO_env.txt` in this folder or the CBIG python 3 environment found under `$CBIG_CODE_DIR/setup/python_env_setup/CBIG_python3_env.yml` (recommended).  

- Preparation prior to AMICO

  Prior to AMICO, please place the diffusion images, bval and bvec files in a single folder. Diffusion files should be named as `${subject_id}.nii.gz`
  Bval and bvec files should have the same root name as the diffusion files (i.e. `${subject_id).bval`). Additionally, if you have a path of b=0 mask
  files, you can specify it in order to skip the b=0 mask calculation step. Each subject should have it's own folder in the mask directory and be named 
  as `${subject_id}_bet_b0_mask.nii.gz`. The file structure for the diffusion image directory should be as follows: 

```
Diffusion input directory (dwi_dir)
|-- sub-1
|-- sub-2
|-- sub-3
    |-- sub-1.nii.gz
    |-- sub-2.bval 
    |-- sub-3.bvec

mask directory (mask_output_dir)
|-- sub-1
|-- sub-2
|-- sub-3
    |-- sub-1_bet_b0_mask.nii.gz

```

- Step 2: Running the script

  Prepare a list of subjects to run. The script can then be run with the following command. Remember to specify the same python environment as 
  your installation:
```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_diffusion_processing2022/AMICO/CBIG_DiffProc_runAMICO.sh \
        --subj_list /path/to/txtfile --dwi_dir /path/to/dwi_images \
        --output_dir /path/to/output --py_env name_of_AMICO_environment \
        --mask_output_dir /path/to/stored/b0_brainmask
```

- Step 2: Results

   An output directory will be generated. You can check the `log` file to see whether there were any errors encountered by each subject. 
   Once completed, the results for each subject can be found under `output/AMICO/NODDI/subject`.


----

## Bugs and Questions

Please contact Leon Ooi at leonooiqr@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.