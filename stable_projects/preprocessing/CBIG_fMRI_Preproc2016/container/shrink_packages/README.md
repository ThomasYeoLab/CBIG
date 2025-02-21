# Shrink Packages

This folder contains scripts to remove redundant files in packages/software used in our fMRI preprocessing pipeline.

## Steps to shrink packages

The original packages (AFNI, ANTs, Freesurfer and FSL) are too large to be entirely copied to the image. We shrink the packages by removing files that we do not use during the preprocessing pipeline.

To get the list of files that are used during the execution of the preprocessing pipeline, we can use the `strace` command. However, since `strace` does not work in the containers, we have to run the unit tests on the HPC with environments/configurations similar to that of the containers (e.g., use the MATLAB Runtime instead of the full MATLAB) and strace the matlab process to get the list of accessed/used files.

As a start, to create an environment similar to containers', we need to:
-   Download and install matlab runtime
-   Copy the original packages to the `../apps` directories
-   Install a conda environment using the `CBIG_python3_container_env.yml` file in the upper directory `../` and activate this environment instead of the default environment in the `~/.bashrc` file
-   Copy/symlink the `./CBIG_preproc_shrink_package_on_HPC_config.sh` to the `~/setup` and source this config instead of the default one in the `~/.bashrc` file

Now that we have setup the environment, we can `strace` the unit test to get the accessed/used files during preprocessing: `cd ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject; source ~/setup/CBIG_preproc_shrink_package_on_HPC_config.sh && strace -ft matlab -nodisplay -nosplash -nodesktop -r "runtests('CBIG_preproc_single_subject_unit_test'); exit;" > unit_test.log 2> strace.log`

Then we can iteratively perform the following steps to shrink the packages:

1.  run one round of unit testing using the above commands to see if the setup works
2.  if the current setup works, use the `./CBIG_preproc_shrink_package.sh`, which relies on the previously generated `strace.log` file, to iteratively move a folder in the `../apps/` directory to a temporary directory such as `../tmp`. If the new setup still works, continue to move the next folder. If the new setup does not work, move the folder back to the `../apps/` directory using the `./CBIG_preproc_shrink_package.sh` script and try another folder.
