# Code for phenotypes and meta-set processing 

----

## References

+ He, T., An, L., Feng, J., Bzdok, D., Eickhoff, S.B. and Yeo, B.T., 2020. [**Meta-matching: a simple approach to leverage large-scale brain imaging datasets to boost prediction of non-imaging phenotypes in small datasets**](https://doi.org/10.1101/2020.08.10.245373), under review.

----

## Usage

This is the code for phenotypes and meta-set processing. In order to use it, you need to have access to UK Biobank and download the csv (https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide) and 25753-2.0 (https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=25753), they correpsond to the `CSV_RAW` and `DIR_PFC` entry in `$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing/config.py`. `CSV_RAW` is the csv directly downloaded from UK Biobank. `DIR_PFC` should contains txt file for each individual participants' 25753-2.0 data field, each txt file should follow xxxxxxx_25753_2_0.txt (where xxxxxxx is eid for this participant). For CBIG internal user, current `$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing/config.py` have the `CSV_RAW` and `DIR_PFC` entry set to files saved in `$CBIG_REPDATA_DIR`, you do not need to change anything. For external user, please update them to the location you downloaded from UK Biobank.

Our code is based on the UK Biobank data field that we have access to. It is very likely you have a different sets of data field that you have access to. And the data field will also be different depends on the time you downloaded it. Therefore, please update the code according to your own data.

### Data for Experiment 1 (UK Biobank only)
* Steps (1 to 4) to run the code to get data for experiment 1
```bash
data_dir="$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing"

cd $data_dir/step1_coarse_filter
# For CBIG user, this line of code need to run under head node
sh CBIG_MM_data_step1_coarse_filter_submit.sh
# first step performs coarse filter to phenotypes of UK Biobank dataset

# wait the step1 code finish
cd $data_dir/step2_krr_filter
# For CBIG user, this line of code need to run under head node
sh CBIG_MM_KRR_filter_wrapper.sh
# this step performs KRR to filter phenotypes with relative good prediction

# wait all jobs finished
python3 CBIG_MM_KRR_filter_summary.py
# summary the KRR result from last step 

cd $data_dir/step3_post_krr_pca
python3 post_krr_pca.py
# this step runs PCA to part of phenotypes filtered to remove redundant ones

cd $data_dir/step2_krr_filter
# For CBIG user, this line of code need to run under head node
sh CBIG_MM_KRR_filter_wrapper.sh pca
# this step performs KRR to filter PCAed phenotypes

# wait all jobs finished
python3 CBIG_MM_KRR_filter_summary.py --pca
# summary the KRR result from last step 

cd $data_dir/step4_experiment_1
python3 experiment_1.py
# this step prepare all the data needed for experiment 1 of the paper
```

### Data for Experiment 2 (UK Biobank and HCP)
Make sure you run the previous steps, here is step to get data for experiment 2.
#### Step 5.1 (step 5 for whole data processing code)
generate subjects/file list for image (seperate for different scan) after ICA-FIX for all subjects
```sh
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing/step5_hcp_data
python3 CBIG_MM_HCP_generate_TC_subj_list.py
```

#### Step 5.2
Extract time series for all subjects from the individual image (seperate for different scan) after ICA-FIX 
```matlab
base_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects/predict_phenotypes/He2022_MM/data_processing/step5_hcp_data');
cd(base_dir)
output_dir = fullfile(base_dir, 'output');
% SUBJ_LIST_HCP_TC in config.py
SUBJ_LIST_HCP_TC = fullfile(output_dir, 'surf_file_list_S1200_1094_210110.txt');
HCP_TC_DIR = fullfile(output_dir, 'TC_419');
ROI_400_NII = fullfile(base_dir, 'extra', 'Schaefer2016_400Parcels_17Networks_colors_19_09_16_subcortical.dlabel.nii');
CBIG_MM_HCP_extract_TC_419(SUBJ_LIST_HCP_TC, HCP_TC_DIR, ROI_400_NII);
```

#### Step 5.3
For each subjects combine the individual time series data, and generate functional connectivity data
```matlab
base_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects/predict_phenotypes/He2022_MM/data_processing/step5_hcp_data');
cd(base_dir)
output_dir = fullfile(base_dir, 'output');
HCP_TC_DIR = fullfile(output_dir, 'TC_419');
HCP_FC_DIR = fullfile(output_dir, 'FC_419');
CBIG_MM_HCP_generate_FC(HCP_TC_DIR, HCP_FC_DIR);
```

#### Step 5.4
generate subjects list for subjects with FC data
```sh
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing/step5_hcp_data
python3 CBIG_MM_HCP_generate_FC_subj_list.py
```

#### Step 5.5
for each subjects from last step, check whether they have all the 58 behaviors that we used. If yes, keep them, if not, remove them. Generate the list of subjects and behaviors after filtering.
```matlab
base_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects/predict_phenotypes/He2022_MM/data_processing/step5_hcp_data');
cd(base_dir)
output_dir = fullfile(base_dir, 'output');
HCP_dir = '/mnt/isilon/CSC1/Yeolab/Data/HCP/S1200/scripts';
measure_list_dir = fullfile(base_dir, 'extra', 'measures_lists');
CBIG_MM_HCP_save_data(output_dir, HCP_dir, measure_list_dir);
```

#### Step 5.6
Run 100 times 10 fold CV with 100 random split. And only keep behaviors with more than 0.1 correlation. Save out the phenotype list and phenotypes.
```matlab
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
base_dir = fullfile(CBIG_CODE_DIR, 'stable_projects/predict_phenotypes/He2022_MM/data_processing/step5_hcp_data');
cd(base_dir)
output_dir = fullfile(base_dir, 'output');
HCP_dir = '/mnt/isilon/CSC1/Yeolab/Data/HCP/S1200/scripts';
measure_list_dir = fullfile(base_dir, 'extra', 'measures_lists');
CBIG_MM_HCP_krr_cv_filter(CBIG_CODE_DIR, HCP_dir, measure_list_dir, output_dir);
```

#### Step 6.0 (optional)
For CBIG internal user, please copy the functional connectivity with 419 ROI for UK biobank data from NSCC `/home/projects/11000481/group/data/uk_biobank/scripts/data_prepare/5_convert_mni2fc/ukbb_fc_419` to `$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing/step6_experiment_2`. It is optional because it takes 20TB more space to store the rsfMRI data and intermediate data to process them.

For external user, you can refer following code to convert UK biobank rsfMRI data data-field 20227 to functional connectivity with 419 ROI. Please make sure you have enough storage. The ~40,000 participants rsfMRI from UK biobank takes more than 20TB to store and process.
* Convert rsfMRI to MNI152 2mm space
```sh
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing/step6_experiment_2/ukbb_20227_to_fc419
in_dir="/ukbb_20227" # folder with downloaded uk biobank 20227 data field, contains file end with '*20227_2_0.zip'
fsl_dir="/apps/fsl/5.0.8" # fsl folder, we tested 5.0.8
current_dir="$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing/step6_experiment_2/ukbb_20227_to_fc419"
out_dir="$current_dir/output/ukbb_20227_mni2mm"
python3 1_conv_fmri_2_mni_parallel.py -i $in_dir -f $fsl_dir -c $current_dir -o $out_dir
```

* Generate mask for Schaefer2018 419 FC in MNI152 2mm
```matlab
base_dir = fullfile(getenv('CBIG_CODE_DIR'),...
    'stable_projects/predict_phenotypes/He2022_MM/data_processing/step6_experiment_2/ukbb_20227_to_fc419');
cd(base_dir)
output_dir = fullfile(base_dir, 'output');
fc400_nii = fullfile(base_dir, 'extra',...
    'Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.nii.gz');
CBIG_MM_create_FC419_MNI2mm(output_dir, fc400_nii);
```

* Generate 419 FC
```sh
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing/step6_experiment_2/ukbb_20227_to_fc419
in_dir="$current_dir/output/ukbb_20227_mni2mm"
out_dir="$current_dir/output/ukbb_fc419"
mask_nii="$current_dir/output/FC419_MNI2mm.nii.gz"
python3 2_conv_mni_2_fc419_parallel.py -i $in_dir -m $mask_nii -o $out_dir
```

#### Step 6
write all the data for experiment 2
```sh
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing/step6_experiment_2
python3 experiment_2.py
```

The final processed data is going to be stored in `data_processed` (you can change it at `DIR_FINAL_OUT` of `config.py`).

----

## Bugs and Questions
Please contact He Tong at hetong1115@gmail.com.