# Replication for comparison of Kernel Regression and Deep Neural Network for RSFC based behavioral prediction (KRDNN)

The replication of KRDNN contains three parts:
1. Kernel regression with cross validation (MATLAB)
2. Kernel regression with normal training, validation and testing (MATLAB)
3. Deep neural network (Python)

For MATLAB code, we have tested in R2014a and R2018b.

----

References
==========
+ He, T., Kong, R., Holmes, A., Nguyen, M., Sabuncu, M., Eickhoff, S.B., Bzdok, D., Feng, J. and Yeo, B.T., 2019. [**Deep Neural Networks and Kernel Regression Achieve Comparable Accuracies for Functional Connectivity Prediction of Behavior and Demographics**](https://www.biorxiv.org/content/10.1101/473603v1), under review.

----

Data
====
The data of part 1 and part 2 are HCP dataset. The data of part 3 are UK Biobank dataset.

----

Run
====

### Kernel regression with cross validation (MATLAB)
----
Update the input variable `unstricted_csv`, `stricted_csv`, `FD_file`, `FC_file` in the following MATLAB code, based on the help text of `/../KR_HCP/CBIG_KRDNN_KRR_HCP.m` function. And run from this `replication` folder.
```MATLAB
clear
addpath(genpath(fullfile(pwd, '/../KR_HCP/')));

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR'); % you may also want to update this if you did not set environment for our CBIG repository
subject_list = fullfile(pwd, 'input', 'kr_hcp', 'He2019_hcp_953_subject_list.txt');
unrestricted_csv = '/share/users/imganalysis/yeolab/data/HCP/S1200/scripts/subject_measures/unrestricted_jingweili_12_7_2017_21_0_16_NEO_A_corrected.csv';
restricted_csv = '/share/users/imganalysis/yeolab/data/HCP/S1200/scripts/restricted_hcp_data/RESTRICTED_jingweili_4_12_2017_1200subjects.csv';
FD_file = '/mnt/eql/yeo1/CBIG_private_data/replication/stable_projects/predict_phenotypes/He2019_KRDNN/HCP/FD_subject_953.txt';
FC_file = '/mnt/eql/yeo1/CBIG_private_data/replication/stable_projects/predict_phenotypes/He2019_KRDNN/HCP/FC_subject_953.mat';
output_dir = fullfile(pwd, 'output_kr_hcp');

CBIG_KRDNN_KRR_HCP(CBIG_CODE_DIR, subject_list, unrestricted_csv, restricted_csv, FD_file, FC_file, output_dir)
rmpath(genpath(fullfile(pwd, '/../KR_HCP/')));
```
The result is saved at `output_kr_hcp/final_result.mat` or the location you stated for `output_dir`.

----

### Kernel regression with normal training, validation and testing (MATLAB)
----
1. Update the setup text file `setup_kr_ukbb.txt` or create a new text file based on your dataset location. You can update the setup file based on the description of `CBIG_KRDNN_KRR_UKBB` function.
2. Once updated or created the setup text file, run following MATLAB code from this `replication` folder:
```MATLAB
addpath(genpath(fullfile(pwd, '/../KR_UKBB/')));
setup_file = fullfile('input', 'kr_ukbb', 'setup_kr_ukbb.txt');
lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20 30 40 50 60 70 80 100 150 200];
CBIG_KRDNN_KRR_UKBB(setup_file, 1000, lambda_set)
rmpath(genpath(fullfile(pwd, '/../KR_UKBB/')));
```
The result is saved at `output_kr_ukbb\final_result.mat` or the location you stated in setup file.

----

### Deep neural network (Python)
----
1. Data preparation
After running previous 2 MATLAB parts, update the first 5 lines in the following sh code and run the following code in terminal from this `replication` folder.
* data_dir: the directory to save the intermediate data for DNN, it could be any where, but it is better to have large space (may need 100GB)
* hcp_fc: path of HCP FC file, which should be same as `FC_file` in part 1
* hcp_kr_output_dir: the output directory for part 1 code. If you did not modify it in the part 1, the default one should work well
* ukbb_fc: path of UKBB FC file, which should be same as line 5 of `setup_kr_ukbb.txt` in part 2
* ukbb_kr_output_dir: the output directory for part 2 code. If you did not modify it in the part 2, the default one should work well
```sh
data_dir=$PWD/data_dnn
hcp_fc=/mnt/eql/yeo1/CBIG_private_data/replication/stable_projects/predict_phenotypes/He2019_KRDNN/HCP/FC_subject_953.mat
hcp_kr_output_dir=$PWD/output_kr_hcp
ukbb_fc=/mnt/eql/yeo1/CBIG_private_data/replication/stable_projects/predict_phenotypes/He2019_KRDNN/UKBB/ukbb_ht_180205_FC_55.mat
ukbb_kr_output_dir=$PWD/output_kr_ukbb
mkdir $data_dir
mkdir $data_dir/original_data_953
mkdir $data_dir/original_data_ukbb_8868
rsync -vP $hcp_fc $data_dir/original_data_953
rsync -vP $PWD/input/kr_hcp/He2019_hcp_953_split.mat $data_dir/original_data_953
rsync -axD $hcp_kr_output_dir/y $data_dir/original_data_953
rsync -vP $ukbb_fc $data_dir/original_data_ukbb_8868
rsync -vP $ukbb_kr_output_dir/ukbb_subject_split.mat $data_dir/original_data_ukbb_8868
rsync -axD $ukbb_kr_output_dir/y $data_dir/original_data_ukbb_8868
```
2. Process data
Update `BASE_DIR` in `/../cbig/He2019/config.py` to location of `data_dir` in previous sh command.
Run the following code in terminal from this `replication` folder:
```sh
export PYTHONPATH=${PYTHONPATH}:$PWD/../
python3 CBIG_KRDNN_proc_data.py
```
3. Run DNN
Run each python3 command in `CBIG_KRDNN_DNN_command.sh` from this `replication` folder (The hyperparameter in each python3 command are final hyperparameters for our results in our paper, tuned by hand-tuning (HCP dataset) and automatic hyperparameter tuning algorithm ([HORD](https://github.com/ilija139/HORD), UK Biobank dataset)). Since it takes a long time and you may need to choose which gpu to run the code on, it is better to run these commands individually. You can use the '--gpu' flag to choose which gpu in your machine to run the code.
4. Result reference
You can find my reference result in `CBIG_KRDNN_DNN_command.sh` for each command
----

Bugs and Questions
====
Please contact He Tong at hetong1115@gmail.com.

