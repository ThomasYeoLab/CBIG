# Code for Meta-matching replication

----

## References
+ He, T., An, L., Chen, P., Chen, J., Feng, J., Bzdok, D., Holmes, A.J., Eickhoff, S.B. and Yeo, B.T., 2022. [**Meta-matching as a simple framework to translate phenotypic predictive models from big to small data**](https://doi.org/10.1038/s41593-022-01059-9), Nature Neuroscience 25, 795-804.

----

## Usage

If you are not an internal CBIG user, you may want to check out the data processing code in `$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing`. After processing, please update the `input_dir` in the `CBIG_MM_KRR_XXX_wrapper_XXX.sh` to your output of data processing code in `$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing`.

### 1. Classical KRR
#### Experiment 1 UK Biobank dataset
Run the classical KRR on UK biobank dataset by:
```bash
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/replication
sh CBIG_MM_KRR_classical_wrapper_UKBB.sh
```

After all jobs finished in previous part, run following code in matlab:
```matlab
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
base_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'He2022_MM');
code_dir = fullfile(base_dir, 'KRR_CLASSICAL');
cd(code_dir)
data_dir = fullfile(base_dir, 'replication', 'output_KRR_classical_ukbb');
phe_list = fullfile(getenv('CBIG_REPDATA_DIR'), 'stable_projects',...
    'predict_phenotypes', 'He2022_MM', 'ukbb_test_final_phe_list.txt');
n_rng = 100;
CBIG_MM_KRR_classical_summary(CBIG_CODE_DIR, data_dir, phe_list, n_rng)
```

#### Experiment 2 HCP dataset
Run the classical KRR on HCP dataset by:
```bash
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/replication
sh CBIG_MM_KRR_classical_wrapper_HCP.sh
```

After all jobs finished in previous part, run following code in matlab:
```matlab
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
base_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'He2022_MM');
code_dir = fullfile(base_dir, 'KRR_CLASSICAL');
cd(code_dir)
data_dir = fullfile(base_dir, 'replication', 'output_KRR_classical_HCP');
phe_list = fullfile(getenv('CBIG_REPDATA_DIR'), 'stable_projects',...
    'predict_phenotypes', 'He2022_MM', 'HCP_diff_roi_final_phe_list.txt');
n_rng = 100;
CBIG_MM_KRR_classical_summary(CBIG_CODE_DIR, data_dir, phe_list, n_rng, 'HCP')
```

### 2. KRR Meta-matching (experiment 1)
Run the KRR meta-matching for experiment 1 by:
```bash
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/replication
sh CBIG_MM_KRR_MM_wrapper.sh
```

After all jobs finished in previous part, run following code in matlab:
```matlab
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
base_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'He2022_MM');
code_dir = fullfile(base_dir, 'KRR_MM');
cd(code_dir)
data_dir = fullfile(base_dir, 'replication', 'output_KRR_mm');
krr_classical_dir = fullfile(base_dir, 'replication', 'output_KRR_classical_ukbb');
input_dir = fullfile(getenv('CBIG_REPDATA_DIR'), 'stable_projects', 'predict_phenotypes', 'He2022_MM');
mm_rng_nums = 100;
CBIG_MM_KRR_MM_summary(CBIG_CODE_DIR, base_dir, data_dir, krr_classical_dir, input_dir, mm_rng_nums)
```

### 3. Run DNN

#### Get spilt files from previous section for DNN
The DNN meta-matching uses the split from KRR classical and meta-matching, if you want to run it, you need to run following command to get the split:
```bash
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/replication

# for experiment 1 UKBB only meta-matching
python3 ../cbig/He2022/get_split.py

# for experiment 2 HCP and UKBB cross datasets meta-matching
python3 ../cbig/He2022/get_split.py -d HCP
```

##### DNN data
Similar as KRR, for CBIG internal user, please copy `ukbb_dnn_input_test.npz` (experiment 1) and `ukbb_dnn_input_cross_dataset.npz` (experiment 2) from `$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/He2022_MM` to `$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/replication/input_dnn` for the input of DNN. 

For external user, please check the data processing code in `$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing`. After running data processing, please copy `NPZ_FNN_INPUT` (experiment 1) and `NPZ_FNN_INPUT_ACROSS_DATASET` (experiment 2) entry based on your location set in `$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing/config.py` to `$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/replication/input_dnn` folder for the input data of DNN.

#### Run DNN
Run the python3 command in `CBIG_MM_DNN_command.sh` from this `replication` folder for experiment 1 and 2 (The hyperparameter in each python3 command are final hyperparameters for our results in our paper, tuned by hand-tuning (HCP dataset) and automatic hyperparameter tuning algorithm ([HORD](https://github.com/ilija139/HORD), UK Biobank dataset)). For the code with "--gpu" flag, it takes a long time and you may need to choose which gpu to run the code on.

You may want to update the directory in `../cbig/He2022/config.py` based on your usage.

### 4. Result reference

#### KRR result reference
* Result for Classical KRR experiment 1 UK Biobank dataset
```matlab
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
base_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'He2022_MM');
result_file = fullfile(base_dir, 'replication', 'output_KRR_classical_ukbb', 'final_result', 'krr_classical_res_test.mat');
result = load(result_file);
disp(squeeze(mean(mean(result.meta_cor, 2), 1)))
% reference disp result:    0.0378    0.0525    0.0759    0.0968    0.1205
% average correlation for K = 10        20        50        100       200
disp(squeeze(mean(mean(result.meta_cod, 2), 1)))
% reference disp result:    -0.0003   -0.0005    0.0021    0.0104    0.0231
% average COD for         K = 10        20        50        100       200
```

* Result for Classical KRR experiment 2 HCP dataset
```matlab
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
base_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'He2022_MM');
result_file = fullfile(base_dir, 'replication', 'output_KRR_classical_HCP', 'final_result', 'krr_classical_res_test.mat');
result = load(result_file);
disp(squeeze(mean(mean(result.meta_cor, 2), 1)))
% reference disp result:    0.0300    0.0475    0.0788    0.1116    0.1562
% average correlation for K = 10        20        50        100       200
disp(squeeze(mean(mean(result.meta_cod, 2), 1)))
% reference disp result:    -0.0135   -0.0156   -0.0075    0.0053    0.0260
% average COD for         K = 10        20        50        100       200
```

* Result for KRR MM experiment 1
```matlab
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
base_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'He2022_MM');
result_file = fullfile(base_dir, 'replication', 'output_KRR_mm', 'meta_result', 'meta_res_test.mat');
result = load(result_file);
disp(squeeze(mean(mean(result.meta_cor, 2), 1)))
% reference disp result:    0.0834    0.1085    0.1357    0.1510    0.1641
% average correlation for K = 10        20        50        100       200
disp(squeeze(mean(mean(result.meta_cod, 2), 1)))
% reference disp result:    -0.0657   -0.0350    0.0045    0.0251    0.0387
% average COD for         K = 10        20        50        100       200
```

#### DNN result reference
You can find my reference result in `CBIG_MM_DNN_command.sh` under each command.

### 5. Haufe transform
After you finished previous steps (espcially the experiment 2 part), you may want to perfrom the Haufe transform.
* First you need to get data for predictive network features (PNF) calculation for "ground truth" and KRR 100 shot by
```bash
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/replication
sh CBIG_MM_KRR_classical_wrapper_haufe_all.sh # "ground truth"

cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/replication
sh CBIG_MM_KRR_classical_wrapper_haufe_100.sh # classical KRR on 100 shot
```

* Then save the result in numpy npz format
```bash
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/replication
python3 ../cbig/He2022/CBIG_haufe_data.py
```

* Run following code for PNF for three different meta-matching case, please make sure you have enough storage space (> 300GB `LARGE_DATA_DIR` in `../cbig/He2022/config.py`) and ram (>= 256GB) to run these.
```bash
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/replication
# PNF: Basic Meta-matching (DNN) Training
python3 ../cbig/He2022/CBIG_ukbb_dnn_haufe.py --gpu 0
# PNF: Basic Meta-matching (DNN)
python3 ../cbig/He2022/CBIG_ukbb_dnn_mm.py --across-dataset --haufe-save
# PNF: Advanced Meta-matching (stacking)
python3 ../cbig/He2022/CBIG_ukbb_dnn_mm_stacking.py --across-dataset --haufe-save
```

* Run the scripts for haufe transform, make sure you have enough storage space (> 300GB LARGE_DATA_DIR in ../cbig/He2022/config.py) and ram (>= 256GB) to run these. You may need to rerun the script multiple times in case it failed due to memory error, the program will save out the intermediate result in case of failing.
```bash
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/replication
python3 ../cbig/He2022/CBIG_haufe_check.py
```

----

## Bugs and Questions
Please contact Tong He at hetong1115@gmail.com, Lijun An at anlijun.cn@gmail.com and Pansheng Chen at chenpansheng@gmail.com.
