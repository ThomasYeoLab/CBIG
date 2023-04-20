References 
=====================
+ Yan, Xiaoxuan, Ru Kong, Aihuiping Xue, Qing Yang, Csaba Orban, Lijun An, Avram J. Holmes et al. [**Homotopic local-global parcellation of the human cerebral cortex from resting-state functional connectivity**](https://doi.org/10.1016/j.neuroimage.2023.120010) NeuroImage (2023): 120010.

Replication for level 400 hMRF parcellation
============================================

We replicate the final 400-level parcellation in our paper using resting-state data of 1479 GSP subjects.
Notice that all filenames and directories in this replication **works for CBIG lab only**.

Input data
==========

**Input**
+ `./input/lh_GSP_subject_fullset_fullpath.csv`: fullpath to all fMRI gifti files of 1479 subjects, left hemisphere only.
+ `./input/rh_GSP_subject_fullset_fullpath.csv`: fullpath to all fMRI gifti files of 1479 subjects, left hemisphere only.

How to run
==========

On the headnode, open a terminal under the current folder. In bash, first run:
`source ./config/CBIG_hMRF_tested_config.sh`

Next, run:
`sh CBIG_hMRF_replicate_400level_parcellation_wrapper.sh ${your_output_dir}`
The above script calls 3 wrapper functions:

- Step 1 calls the script `CBIG_hMRF_generate_subject_fullpath_GSP.sh` and generate the full subject lists.
- Step 2 calls the code under `../code/step1_generate_fmri_input` and generate the premultiplied fMRI matrix. 
- Step 3 calls the code under `../code/step2_generate_parcellation` and generate the level-400 parcellation.

The total time for replicating the 400-level parcellation should be within 10 hours.
To check the replication results, open Matlab and `cd` to the current directory, run `CBIG_hMRF_check_replication_results(your_output_dir)`.

Replicating results for other resolutions (100 to 1000, excluding 400)
============================================================
Please refer to the wrapper files under `${CBIG_REPDATA_DIR}/stable_projects/brain_parcellation/Yan2023_homotopic/all_resolution_wrapper` to replicate the parcellations for all resolutions other than the 400 version.