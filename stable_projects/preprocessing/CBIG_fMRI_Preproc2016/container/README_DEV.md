# Containerized CBIG_fMRI_Preproc2016 Developer Guide

This folder contains the scripts to containerize the CBIG fMRI 2016 preprocessing pipeline.

## Build Docker Image

### Shrink packages

Since the original packages (AFNI, ANTs, Freesurfer and FSL) are too large to be fully incorporated to the images. We have developed a script to remove files that are redundant to our preprocessing pipeline. The detailed description of the steps to shrink the packages can be found in the [shrink_packages README](shrink_packages/README.md).

### Run `docker build`

To build docker image:

```
docker build -t thomasyeolab/cbig_fmri_preproc2016:latest \
-f $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/container/Dockerfile \
$CBIG_CODE_DIR --progress=plain 2>&1 | tee build.log
```

Here, a `build.log` file will be created to store the docker build logs. Moreover, the file after the `-t` flag is the name of the image, and the file after the `-f` flag is the path to the Dockerfile, which contains the instructions to build the image.

### Changes made to the existing files

During the build process as specified in the Dockerfile, other than installing the required packages and setting up the environments, all CBIG_fMRI_preproc2016 related files will be copied to the image, with several modifications to the following files:

-   `CBIG_fMRI_Preproc2016/CBIG_preproc_fMRI_preprocess.csh`: changed the default matlab to use the matlab runtime by `set matlab_runtime_util = "$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/container/utilities"`. And removed the -z flag from rsync commands as compression may result in errors in singularity containers.
-   `CBIG_pbsubmit`: changed the command to directly run the command instead of submitting to the job scheduler
-   `CBIG_sample_config.sh`: used the correct path that points to the `CBIG_private` repository. Properly configure conda, activate CBIG_py3, and export variables required by the preprocessing pipeline
-   `CBIG_generic_setup.sh`: commented out the symlinking pre-comit hook and do not setup Workbench directory
-   `CBIG_preproc_single_subject_unit_test.m`: remove lines that ssh into headnode when performing the unit tests

## Unit Testing

To run unit test on the containerized images, we have to run it on the HPC since the licensed MATLAB and the test data are stored on it.

-   To run the unit test using a docker container:

```

export MATLAB_RUNTIME_DIR=/apps/matlab/MATLAB_Compiler_Runtime/v95 && docker run -it --rm -v $CBIG_MATLAB_DIR:$CBIG_MATLAB_DIR:ro -v $CBIG_TESTDATA_DIR:$CBIG_TESTDATA_DIR:ro -v $MATLAB_RUNTIME_DIR:/apps/matlab/MCR_R2018b/v95:ro -v ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/output:/app/CBIG_private/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/output --entrypoint /bin/bash thomasyeolab/cbig_fmri_preproc2016:latest -c "source ~/.bashrc && cd \$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject && matlab -nodisplay -nosplash -nodesktop -r \"runtests('CBIG_preproc_single_subject_unit_test'); exit;\""

```

-   To run the unit test using a singularity image, first load the singularity module by `module load singularity/3.11.1`, then run the following command:

```

export MATLAB_RUNTIME_DIR=/apps/matlab/MATLAB_Compiler_Runtime/v95 && singularity exec --cleanenv --no-home --writable -B $CBIG_MATLAB_DIR:$CBIG_MATLAB_DIR:ro -B $CBIG_TESTDATA_DIR:$CBIG_TESTDATA_DIR:ro -B $MATLAB_RUNTIME_DIR:/apps/matlab/MCR_R2018b/v95:ro -B ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/output3:/app/CBIG_private/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/output cbig_fmri_preproc2016 bash -c "source /app/setup/CBIG_config_in_container.sh && cd /app/CBIG_private/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject && matlab -nodisplay -nosplash -nodesktop -r \"runtests('CBIG_preproc_single_subject_unit_test'); exit;\""

```

In addition to the above, to truly test the outputs users will get, a unit test is created to test both Singularity and Docker images from a user's perspective. Once the docker image is uploaded to the Docker Hub, you may run the unit test located at `${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_containerization_unit_test.sh`.

## Push to Docker Hub

-   First, generate a [Personal Access Token (PAT)](https://docs.docker.com/security/for-developers/access-tokens/) for the `thomasyeolab` Docker Hub account
-   Then type in `docker login -u thomasyeolab` and paste the generated PAT
-   Finally, push the image to the Docker Hub by running `docker push thomasyeolab/cbig_fmri_preproc2016:latest`

## Setup GitHub Workflow to Sync User Guide README with Docker Hub Description

-   The workflow is already written at `${CBIG_CODE_DIR}.github/workflows/sync_dockerhub_description.yml`, but the Docker Hub PAT (Personal Access Token) is required to be stored in the GitHub Action secrets of the upstream CBIG_private repository
-   Go to Setting tab of upstream CBIG_private repository and click on the "Secrets and variables" tab on the left.
-   Then click on the revealed "Action" tab and click on the "New repository secret" button on the right.
-   Fill in two secrets: `DOCKERHUB_USERNAME` and `DOCKERHUB_PAT` with the Docker Hub username (i.e., thomasyeolab) and the PAT respectively.
-   Then the workflow should automatically run whenever there are changes to files in the `stable_projects/preprocessing/CBIG_fMRI_Preproc2016/container` directory.
