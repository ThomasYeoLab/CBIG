Diffusion processing project depends on real diffusion MRI data, which may vary from case to case. We provide examples on how to run the data for the following:

1. `diffusionQC`: Run `CBIG_DiffProc_diffusionQC_run_example.m` to generate the example data. We provide the input and output for your reference. Please refer to the top level `CBIG2022_DiffProc` folder for specific instructions relating to this script. Data is from group 0 and group 1 of the NKI diffusion data, which can be accessed [here](http://fcon_1000.projects.nitrc.org/indi/enhanced/data/download_nii8.html). Use `CBIG_DiffProc_diffusionQC_check_example_results.m` to check if your results are consistent with the example.

2. `TBSS`: Run `CBIG_DiffProc_TBSS_run_example.m` to generate the example data. We provide the input and output for your reference. Please refer to the `TBSS` folder for specific instructions relating to this script. Data is from the resulting FA and MD maps generated from `CBIG_DiffProc_diffusionQC_run_example.m`. Use `CBIG_DiffProc_TBSS_check_example_results.m` to check if your results are consistent with the example. Additional Note: The script has been run without changing the default thresholds - the TBSS skeleton could be further improved by tuning the registration and the FA thresholds

3. `AMICO`: Run `CBIG_DiffProc_AMICO_run_example.m` to generate the example data. We provide the input and output for your reference. Please refer to the `AMICO` folder for specific instructions relating to this script. Data is from the NODDI example dataset found [here](https://github.com/daducci/AMICO/wiki/NODDI). Use `CBIG_DiffProc_AMICO_check_example_results.m` to check if your results are consistent with the example.

Please note the following: 

a) There is no example for MRtrix processing as the resulting files are too big to be stored on Github.

b) The examples are different from the data used in the `unit tests`, as the data in the unit tests allows testing for more comprehensive scenarios, however, they are not publicly available. Therefore, the `unit_tests` example can only run within CBIG lab.
