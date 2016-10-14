Please read the following instructions **carefully** in order to have your local environment compatible with CBIG repository.

# 1) Copy and modify configuration file
`CBIG` repository uses many external softwares. In order to easily use and manage these softwares, we use a script to configure the necessary paths and variables.

A sample configuration script is saved at:
```
<your-cbig-directory>/setup/CBIG_sample_config.sh (for Bash)
<your-cbig-directory>/setup/CBIG_sample_config.csh (for C-shell)
```

Looking into these sample configuration scripts, you can find environmental variables that store paths to directories of various softwares. Note that we point such environmental variables to specific version of each software. We **highly recommend** this practice since different versions of a software can greatly affect outputs of your code.

###a) Copy configuration script

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

###b) Modify your configuration script

Open your newly copied configuration script with a text editor and change the environmental variables to point to the correct directories of the respective softwares. For example: change `CBIG_MATLAB_DIR` to point to your installation of Matlab, and `CBIG_FSLDIR` to your installation of `FreeSurfer`.

###c) Source your configuration script from `.bashrc` or `.cshrc`

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

###Check if your configuration works
- If you are still in BASH or C-SHELL from the previous step, exit.
- Open a new BASH or C-SHELL.
- If you have `FreeSurfer` path configured properly, you should see a few lines from FreeSurfer when your command-line interface is booting up. One of those lines is:
```Setting up environment for FreeSurfer/FS-FAST (and FSL)```
- If you do
```
cd $CBIG_CODE_DIR
```
You should be in the directory of CBIG repository.

# 2) Set `MATLABPATH` and `startup.m`
`MATLABPATH` is a variable containing paths that Matlab will search (when it is launched) to locate files being used by Matlab. 

`startup.m` is a Matlab script that executes a set of commands of your choosing to set up your workspace in Matlab (when it is launched). We use `startup.m` to define paths to functions inside CBIG repository so that you can call them directly in Matlab.

We will add a `startup.m` of CBIG repository to `MATLABPATH` so that whenever you start Matlab, it will look into `MATLABPATH`, locate `startup.m` and set up your Matlab workspace to work with Matlab code under `CBIG` repository.

###a) Set `MATLABPATH` in your `.bashrc` or `.cshrc`
In the following example, we will keep `startup.m` under a folder `matlab` in your home directory.

**[IMPORTANT]**: you may already a have `startup.m` under `~/matlab` folder. If so, please **check carefully** that you are not overwriting important configurations of your existing `startup.m` with our version. We **highly recommend** that you keep an old copy of your `startup.m`.

```
mkdir ~/matlab
```

If you use BASH, add the following line to `.bashrc`:
```
export MATLABPATH=~/matlab
```
If you use C-SHELL, add the following line to `.cshrc`:
```
setenv MATLABPATH=~/matlab
```
###b) Copy `startup.m` from CBIG repository
```
cp $CBIG_CODE_DIR/setup/startup.m ~/matlab
```

