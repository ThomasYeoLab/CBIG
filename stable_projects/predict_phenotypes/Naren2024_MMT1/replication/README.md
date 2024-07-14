# Replication of Meta-matching T1

## References
+ Wulan, Naren, et al. "Translating phenotypic prediction models from big to small anatomical MRI data
  using meta-matching." bioRxiv (2024): 2023-12.
  
----

## Data
If you do not have access to UK Biobank, HCP-YA, HCP-Aging data, you need to apply for these datasets to replicate 
the results in Wulan et al., 2024.

----

## Replication

### 1. Setup

Please setup the conda environment before running replication code.

```
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1

conda env create -f replication/config/CBIG_MMT1_python_env.yml
```

### 2. Data split

Get spilt files for k-shot procedure on meta-test set:

```

ssh headnode

cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1/replication

sh CBIG_get_split.sh 

```
After finishing above experiments, run following command to run DNN models.

### 3. Training DNN

Run following command to train 3D CNN model using meta-training set:

```
ssh headnode

cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1/replication

sh CBIG_train_ukbb_model.sh 

```

In side this file contains training a model for experiment 1 within UKBB 
and a model for experiment 2 HCPYA and experiment 3 HCP-Aging.
After finishing above experiments, run following command to get outputs 
from above trained 3D CNN model using meta-test set:

```
ssh headnode

cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1/replication

sh CBIG_get_output.sh 

```

In side this file contains getting outputs from the pretrained model, the outputs will then be used for classical
transfer learning, meta-matching finetune and meta-matching stacking. 
After finishing above experiments, run following command to perfrom k-shot procedure.

### 4. K-shot procedure

To perform k-shot procedure on meta-test set, you need to run following command:

```

ssh headnode

cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1/replication

# for experiment 1 UKBB only 
sh CBIG_withinUKBiobank_exp.sh 

# for experiment 2 UKBB -> HCPYA 
sh CBIG_HCPYA_exp.sh 

# for experiment 3 UKBB -> HCP-Aging
sh CBIG_HCPA_exp.sh 

```

### 5. Compare with reference results

Once you have finished the replication above, you could run the following command to check whether 
you have successfully replicated results.

```

ssh headnode

cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1/replication

sh CBIG_check_replication_res.sh

```

### 6. Haufe transformation

To perform Haufe transformation, you could run the following command:

```

ssh headnode

cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1/replication

sh CBIG_haufe_check.sh

```

----
### Clean up

Run the following commands to clean up after running experiments

```
rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1/replication/results

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1/job_logs
```

----

## Bugs and Questions
Please contact Naren Wulan at wulannarenzhao@gmail.com, Chen Zhang at chenzhangsutd@gmail.com and 
Lijun An at anlijuncn@gmail.com
