## Replication overview:
* **Step 0.** Initialize project
* **Step 1.** Initialize pipeline
* **Step 2.** Split rs-fMRI into 1-min windows
* **Step 3.** Create simulated runs (by mix-&-match windows)
* **Step 4.** Run Kong2026 tree TMS targeting
  * **Sub-step 4A.** Ran on non-ground-truth (simulated) data
  * **Sub-step 4B.** Ran on ground-truth data
  * Note within tree targeting algorithm, MS-HBM parcellation generated.
* **Step 5.** Compare similarity to ground-truth
* **Step 6.** Compare results to reference outputs

Note that all filenames and directories in this replication **work for CBIG lab only**.

---

## References
Tan, T. W. K., Kong, R., Xue, A., Cheng, J., et al. (2026). [Excessive Censoring Degrades Individual-Specific Cortical Parcellations and Personalized TMS Targets]. (https://doi.org/10.64898/2026.03.09.710457). _bioRxiv_, 2026.2003.2009.710457.

---

## Folder structure
* `config`: Contains conda env .yml, and tested config.
* `init`: Contains Python scripts to initialize project & pipeline.
* `select_subID_window`: Contains Python wrapper to create 1-min windows.
* `select_subID_sim`: Contains Python wrapper to create simulated runs.
* `TMS_Kong2025`: Contains Python wrappers to run TMS on ground-truth & non-ground-truth/simulated data.
* `compare_Kong2025`: Contains Python wrapper to compare similarity to ground-truth.

---

## Data
In this replication, we run analyses on 7 datasets: HBNSSI, IPCAS6, Kirby, MSC, Myconnectome, Newbold, and NSD.
If you do not have access to [NSD](https://naturalscenesdataset.org/) dataset, please apply for access [here](https://cvnlab.slite.page/p/dC~rBTjqjb/How-to-get-the-data).

---

## Environment setup
**Important!** Please follow all steps under "Setup" section from `${TAN2026_TMSMOTION_DIR}/README.md` before replication.

**Please modify `~/.bashrc` file**.
1. Remove/comment any lines which `source` a CBIG config file, to prevent conflicts. For example, remove/comment `source $CBIG_CODE_DIR/setup/CBIG_sample_config.sh`.
2. Add the following line, while replacing `$CBIG_CODE_DIR` with full path: `source $CBIG_CODE_DIR/stable_projects/preprocessing/Tan2026_TMSmotion/replication/config/CBIG_TMSmotion_tested_config.sh`
3. To initialize environment variables properly, please restart terminal, and restart VNC Viewer port.
---

## Run

### Important notes before running
* Here, we replicate both FD pipeline (`pl1`) and FDRMS pipeline (`pl2`).
* Full replication for 1 pipeline takes ~ 2 weeks to run, and requires ~2 TB storage.
  * Hence, replication has 2 modes: `full` and `lite`.
    * `full` mode runs all participants, simulated sessions, and thresholds.
    * `lite` mode runs complete analysis like `full`, but for 2 participants, 2 simulated sessions, and 2 thresholds.
      * `lite` mode will take ~1 hr 30 min to run, and requires ~3 GB storage.
* Replication wrapper runs all steps (including check against reference outputs) in "Replication overview" section above.

### Steps to run
* Below steps for `lite` mode. However, if want `full` mode, please change `lite` to `full` in wrapper command argument.

#### To run FD pipeline (`pl1`)
Open terminal, enter headnode:
```
ssh headnode
```

On same terminal, call replication wrapper:
```
source CBIG_init_conda; conda activate CBIG_TMSmotion;
$CBIG_CODE_DIR/stable_projects/preprocessing/Tan2026_TMSmotion/replication/CBIG_TMSmotion_replication_wrapper.sh lite pl1
```

#### To run FDRMS pipeline (`pl2`)

Open terminal, enter headnode:
```
ssh headnode
```

On same terminal, call replication wrapper:
```
source CBIG_init_conda; conda activate CBIG_TMSmotion;
$CBIG_CODE_DIR/stable_projects/preprocessing/Tan2026_TMSmotion/replication/CBIG_TMSmotion_replication_wrapper.sh lite pl2
```

* If replication successful, will see `Replication successful for pl1/2 (lite)` message on terminal. Else, will see `Replication unsuccessful for pl1/2 (lite)` message.
* If step fails, can find more detailed log file of step at `$CBIG_CODE_DIR/stable_projects/preprocessing/Tan2026_TMSmotion/replication/logs`.

---

## Additional step-specific info
* If ran replication wrapper above, no need run commands below.
* For each command below, please be in headnode terminal to submit jobs (via `ssh headnode`).
* Below commands use `pl1` and `lite` mode as example. May switch args to `pl2` and/or `full` mode if preferred.

### **Step 0.** Initialize project
* **Purpose:** Create `main_config.json`; contains parameters shared across pipelines.
* **Key output:** `$TAN2026_TMSMOTION_DIR/replication/config/main_config.json`.
* **Expected run time:** ~1 min
* **To run:**
```
source CBIG_init_conda; conda activate CBIG_TMSmotion;
python $TAN2026_TMSMOTION_DIR/replication/init/script/CBIG_TMSmotion_init_project.py
```

### **Step 1.** Initialize pipeline
* **Purpose:** Create `pipeline.json`; contains parameters specific to pipeline.
* **Key output:** `$TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json`.
* **Expected run time:** ~1 min
* **To run:**
```
source CBIG_init_conda; conda activate CBIG_TMSmotion;
python $TAN2026_TMSMOTION_DIR/replication/init/script/CBIG_TMSmotion_init_pipeline.py $TAN2026_TMSMOTION_DIR/replication/init/input/pl1.ini lite
```

### **Step 2.** Split rs-fMRI into 1-min windows
* **Purpose:** Prepare 1-min windows to mix-&-match for later step.
* **Key output:** `$TAN2026_TMSMOTION_DIR/replication/pl1/output/select_subID_window/window_list_1.txt`. Each row contains data of one 1-min rs-fMRI window.
* **Expected run time:** ~3 min
* **To run:**
```
source CBIG_init_conda; conda activate CBIG_TMSmotion;
python $TAN2026_TMSMOTION_DIR/replication/select_subID_window/script/CBIG_TMSmotion_select_subID_window_wrapper.py $TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json
```

### **Step 3.** Create simulated runs (by mix-&-match windows)
* **Purpose:** Mix-&-match 1-min windows to generate simulated runs.
* **Key output:** `$TAN2026_TMSMOTION_DIR/replication/pl1/output/select_subID_sim/result/summary/sim_run_list_3_grid.txt`. Each row contains data of 1 censoring threshold from 1 simulated session from 1 participant from 1 motion bin. 
* **Expected run time:** ~30 min (`lite`), ~3 hrs (`full`)
* **To run:**
```
source CBIG_init_conda; conda activate CBIG_TMSmotion;
python $TAN2026_TMSMOTION_DIR/replication/select_subID_sim/script/CBIG_TMSmotion_select_subID_sim_wrapper.py $TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json
```

### **Step 4.** Run Kong2026 tree TMS targeting
* **Purpose:** Runs Kong2026 TMS algorithm on both ground-truth and non-ground-truth/simulated data in parallel.
* **Key output (ground-truth):** `$TAN2026_TMSMOTION_DIR/replication/pl1/output/TMS_Kong2025_dep/<dataset>/<subID>/gt/TMS/target_results/group_weight/<subID>/MSHBM_scalp_target_results/optimal_target_results.mat`.
* **Key output (non-ground-truth/simulated):** `$TAN2026_TMSMOTION_DIR/replication/pl1/output/TMS_Kong2025_dep/<dataset>/<subID>/ngt/<rm_run_th>/<combo>/<sim_sesID>/FD<vary_FD_th>_<vary_DVARS_th>/TMS/target_results/group_weight/<sim_sesID>/MSHBM_scalp_target_results/optimal_target_results.mat`.

* **To run (ground-truth):**
```
source CBIG_init_conda; conda activate CBIG_TMSmotion;
python $TAN2026_TMSMOTION_DIR/replication/TMS_Kong2025/script/CBIG_TMSmotion_Kong2025_gt_wrapper.py $TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json
```
* **To run (non-ground-truth/simulated):**
```
source CBIG_init_conda; conda activate CBIG_TMSmotion;
python $TAN2026_TMSMOTION_DIR/replication/TMS_Kong2025/script/CBIG_TMSmotion_Kong2025_ngt_wrapper.py $TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json
```

* **Total expected run time (for both ground-truth & non-ground-truth/simulated):** ~40 min (`lite`), ~2 weeks (`full`)


### **Step 5.** Compare similarity to ground-truth
* **Purpose:** For TMS, compute Euclidean distance between simulated target & participant's ground-truth target. For parcellation, compute Dice overlap between simulated parcellation & participant's ground-truth parcellation.
* **Key output:** `$TAN2026_TMSMOTION_DIR/replication/pl1/output/compare_Kong2025_dep/summary/compare_summary_list_4.txt`. Each row contains the Euclidean distance (`eucli_dist`) and Dice overlap (`17_network_dice`) of 1 censoring threshold from 1 simulated session from 1 participant from 1 motion bin.
* **Expected run time:** ~3 min (`lite`), ~5 hrs (`full`)
* **To run:**
```
source CBIG_init_conda; conda activate CBIG_TMSmotion;
python $TAN2026_TMSMOTION_DIR/replication/compare_Kong2025/script/CBIG_TMSmotion_compare_Kong2025_wrapper.py $TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json
```

### **Step 6.** Compare results to reference outputs
* **Purpose:** Compare Euclidean distances (`eucli_dist`) and Dice overlaps (`17_network_dice`) to reference.
* **Expected run time:** ~1 min (`lite`), ~1 hr (`full`)
* **To run:**
```
source CBIG_init_conda; conda activate CBIG_TMSmotion;
python $TAN2026_TMSMOTION_DIR/replication/CBIG_TMSmotion_check_reference_output.py pl1 lite
```

---

## Bugs and questions
Please contact Trevor Tan at [trevortanweikiat97@gmail.com](mailto:trevortanweikiat97@gmail.com) and Thomas Yeo at [yeoyeo02@gmail.com](mailto:yeoyeo02@gmail.com).