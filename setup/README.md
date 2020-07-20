Please read the following instructions **carefully** in order to have your local environment compatible with CBIG repository. 

# 0) Prerequisites
Linux system needs to be installed.
The Linux distribution needs to be supported by the software used by this repository, which are
[FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall),
[FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation),
[SPM](https://www.fil.ion.ucl.ac.uk/spm/),
[AFNI](https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/background_install/install_instructs/index.html),
[ANTs](http://stnava.github.io/ANTs/),
[Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench),and Matlab.
Linux distribution choices include but are not limited to RedHat, SuSE, Ubuntu and their derivatives.

If you want to use a specific stable project, you may not need to install all the above mentioned softwares and may need to download some other packages. Please refer to the project's README page for more information on software requirements.

# 1) Copy and modify configuration file
`CBIG` repository uses many external softwares. In order to easily use and manage these softwares, we use a script to configure the necessary paths and variables.

A sample configuration script is saved at:
```
<your-cbig-directory>/setup/CBIG_sample_config.sh (for Bash)
<your-cbig-directory>/setup/CBIG_sample_config.csh (for C-shell)
```

Looking into these sample configuration scripts, you can find environmental variables that store paths to directories of various softwares. Note that we point such environmental variables to specific version of each software. We **highly recommend** this practice since different versions of a software can greatly affect outputs of your code.

### a) Copy configuration script

Create a new folder to keep your configuration script, for example:
```
mkdir ~/setup
```
- Copy `CBIG_sample_config.sh` if you use `BASH`, or `CBIG_sample_config.csh` if you use `C-SHELL` to your newly created folder:
```
cp <your-cbig-directory>/setup/CBIG_sample_config.sh ~/setup/CBIG_FS5.3_config.sh` (BASH)
```
or
```
cp <your-cbig-directory>/setup/CBIG_sample_config.csh ~/setup/CBIG_FS5.3_config.csh` (C-SHELL)
```
Note that we renamed the config script to clearly reflect the software version we are using. This is useful since we may want to run CBIG repository with different software versions.

**[IMPORTANT]** The above mentioned sample configuration script is a general setup for the whole CBIG repository. Since different projects may have different software configurations, each stable project has its own config file saved under its `config` folder.

For example, a tested configuration script for `Kong2019_MSHBM` projecs is saved at:
```
<your-cbig-directory>/stable_projects/brain_parcellation/Kong2019_MSHBM/config/CBIG_MSHBM_tested_config.sh
```

If you want to use a specific stable project, you need to copy that project's tested config file instead of the CBIG sample config file into your `~/setup` folder 

### b) Modify your configuration script

Open your newly copied configuration script with a text editor and change the environmental variables to point to the correct directories of the respective softwares. For example: change `CBIG_MATLAB_DIR` to point to your installation of Matlab, and `CBIG_FSLDIR` to your installation of `FreeSurfer`.

### c) Source your configuration script from `.bashrc` or `.cshrc`

This step is to run your configuration script every time your shell is launched.

Following the examples in the previous steps, in Bash, add the following line to `~/.bashrc`:
```
source ~/setup/CBIG_FS5.3_config.sh
```
or in C-shell, add the following line to `~/.cshrc`:
```
source ~/setup/CBIG_FS5.3_config.csh
```

**[IMPORTANT]** Please check that in your `.bashrc` or `.cshrc`, there is no sourcing of softwares that conflict with our configuration script.

**[IMPORTANT]** If you use CBIG Python setup, this line must be the last line in your `.bashrc` or `.cshrc`

# 2) Set `MATLABPATH` and `startup.m` 
`MATLABPATH` is a variable containing paths that Matlab will search (when it is launched) to locate files being used by Matlab. 

`startup.m` is a Matlab script that executes a set of commands of your choosing to set up your workspace in Matlab (when it is launched). We use `startup.m` to define paths to functions inside CBIG repository so that you can call them directly in Matlab.

We will add a `startup.m` of CBIG repository to `MATLABPATH` so that whenever you start Matlab, it will look into `MATLABPATH`, locate `startup.m` and set up your Matlab workspace to work with Matlab code under `CBIG` repository.

### a) Set `MATLABPATH` in config file

**[IMPORTANT]**: If you already setup the `MATLABPATH` in your `.bashrc` or `.cshrc`, please **remove** it.

The `MATLABPATH` has already been set in the configuration file as following:
In `CBIG_FS5.3_config.sh`,
```bash
export MATLABPATH=$CBIG_CODE_DIR/setup
```
or
In `CBIG_FS5.3_config.csh`,
```csh
setenv MATLABPATH = $CBIG_CODE_DIR/setup
``` 

### b) Set `startup.m` in CBIG repository
**[IMPORTANT]**: If you have `startup.m` under `~/matlab` folder, please **remove** it.

In our CBIG lab, we use the same `startup.m` in our repo for everyone. This is useful because others can replicate our work easily.
If you want to use some packages written by people outside our lab, you can follow the following steps.
1. Check whether it exists in `external_packages/matlab/default_packages`
If it exists, then you can just use it.
2. Check whether it exists in `external_packages/matlab/non_default_packges`
If it exists, you can use the package by `addpath` of this package in your function, but remember to `rmpath` of this package at the bottom of your function. For the format, see [Convention](https://github.com/YeoPrivateLab/CBIG_private/wiki/Convention).
3. If the package doesn't exist, you can discuss with the admin [minh, ruby, jingwei, nanbo] and make a pull request to add this package into our CBIG repo. See [Level 2 User](https://github.com/YeoPrivateLab/CBIG_private/wiki/Level-2-User:--Contribute-to-CBIG_private-repository) about how to make a pull request. 

# 3) Check if your configuration works
### a) Check FreeSurfer and FSL
- If you are still in BASH or C-SHELL from the previous step, exit.
- Open a new BASH or C-SHELL.
- If you have `FreeSurfer` path configured properly, you should see a few lines from FreeSurfer when your command-line interface is booting up. One of those lines is:
  ```Setting up environment for FreeSurfer/FS-FAST (and FSL)```
- Use `freeview` to load MNI template

```bash
freeview $CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS4.5.0/mri/norm.nii.gz
```

* Use `fsleyes` (for FSL verision 5.0.10 or above) to load MNI template

```bash
fsleyes $CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS4.5.0/mri/norm.nii.gz
```

If you are using FSL with version 5.0.8 or below, you can use `fslview` to load MNI template

```bash
fslview $CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS4.5.0/mri/norm.nii.gz
```

### b) Check Matlab functions
- Open Matlab
- Read a sample volume in MNI space, convert it to fsaverage space and visualize it

```matlab
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
x = MRIread([CBIG_CODE_DIR '/utilities/matlab/figure_utilities/draw_surface_data_as_annotation/sample/sample_vol.nii.gz']);
[lh_data, rh_data] = CBIG_ProjectMNI2fsaverage_Ants(x, 'fsaverage');
CBIG_DrawSurfaceMaps(lh_data, rh_data, 'fsaverage', 'inflated', 1e-5, 5e-5);
```
- You should get similar pattern as  
`$CBIG_CODE_DIR/utilities/matlab/figure_utilities/draw_surface_data_as_annotation/sample/fsaverage_visualization.png`
