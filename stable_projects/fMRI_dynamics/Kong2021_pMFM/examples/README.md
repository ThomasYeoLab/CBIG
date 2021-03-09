# Examples of Parametric Mean Field Model (pMFM) Parameter Estimation
* In this example, conduct the estimation process of pMFM. However, the whole process will take around 30 hours to finish (exact time will depend on GPU model). Therefore, instead of performing the whole process with 500 estimation iterations, the example only performs 5 iterations
* The estimated time for the example is around 20 mins (The GPU model we use for estimation is GTX 1080 Ti)

# Usage
1. In the terminal, go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/examples/scripts`
2. In the terminal, source into the Anaconda environment: `source activate pMFM`
3. In the terminal, run `python CBIG_pMFM_parameter_estimation_example.py`
4. The output of the toy example will be automatically stored into `$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/examples/output` folder named by `example_output.csv`. The folder also contains the expected output from toy example named `expected_output.csv`.
5. To get the same result as `expected_output.csv`, we highly recommend the users to follow the setup instructions to get the proper runtime CUDA version and conda environment.