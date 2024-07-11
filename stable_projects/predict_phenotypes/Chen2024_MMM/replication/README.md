# Code for Multilayer Meta-matching replication

----

## References
+ Chen, P., An, L., Wulan, N., Zhang, C., Zhang, S., Ooi, L. Q. R., ... & Yeo, B. T. (2023). [**Multilayer meta-matching: translating phenotypic prediction models from multiple datasets to small data**](https://www.biorxiv.org/content/10.1101/2023.12.05.569848v1.abstract). bioRxiv, 2023-12.

----
## Data
If you do not have access to [UK Biobank](https://www.ukbiobank.ac.uk/), [ABCD](https://nda.nih.gov/study.html?id=824), [GSP](http://neuroinformatics.harvard.edu/gsp/), [HBN](https://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network), [eNKI](http://fcon_1000.projects.nitrc.org/indi/enhanced/), and [HCP](https://www.humanconnectome.org/) data, you need to apply for these datasets to replicate the results in [Chen et al., 2023](https://www.biorxiv.org/content/10.1101/2023.12.05.569848v1.abstract). 

## Usage
### 0. Set up
Please setup the conda environment before running replication code.
```bash
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM
conda env create -f CBIG_MMM_python_env.yml
source activate CBIG_Chen2024;
```

### 1. Run replication
The wrapper script `replciations/CBIG_MMM_replication.sh` contains all steps to replicate the results in [Chen et al., 2023](https://www.biorxiv.org/content/10.1101/2023.12.05.569848v1.abstract). . For CBIG users, you just need to run the following commands.

```
ssh headnode
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/replication
sh ./CBIG_MMM_replication.sh
```

Note that the hyperparameters used for DNN were automaticaly tuned by [optuna](https://optuna.org/) , you can replicate by runing:

```bash
# run optuna to automatically tune hyper-parameters for DNN training
sh ../scripts/CBIG_MMM_optuna_submit_job.sh
# Accuracy: 0.14080379011481273
# Best hyperparameters: {'dropout': 0.4, 'n_l1': 512, 'n_l2': 256, 'n_l3': 128, 'n_l4': 1024, 'n_layer': 2, 'batch_size': 128, 'lr': 0.036612177992895435, 'weight_decay': 2.727460424379488e-07}
```

### 2. Compare with reference results
Once you have finished the `CBIG_MMM_replication.sh`, you could run the following command to check whether you have successfully replicated results.

```
python ../cbig/Chen2024/check_reference_results.py
```

### 3. Clean up
You may want to run the following commands to clean up all generated files.

```
rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/replication/output

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/replication/output_intermediate

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/replication/log

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/replication/models

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/replication/output_KRR_classical_*
```

----
## Bugs and Questions
Please contact Pansheng Chen at chenpansheng@gmail.com, Lijun An at anlijun.cn@gmail.com, and Chen Zhang at chenzhangsutd@gmail.com.