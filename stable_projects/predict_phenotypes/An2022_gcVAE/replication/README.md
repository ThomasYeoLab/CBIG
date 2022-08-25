# Replication of Goal-specific brain MRI harmonization

## References
+ An, L., Chen, J., Chen, P., Zhang, C., He, T., Chen, C., Zhou, J., Yeo, B.T., 2022. [Goal-specific brain MRI harmonization](https://doi.org/10.1016/j.neuroimage.2022.119570), NeuroImage, In press

----

## Data

If you do not have access to [ADNI](http://adni.loni.usc.edu/), [AIBL](https://aibl.csiro.au/) and [MACC](http://www.macc.sg/) data, you need to apply for these datasets to replicate the results in [An et al., 2022](https://doi.org/10.1016/j.neuroimage.2022.119570). 

For ADNI data, we used the spreadsheets by [TADPOLE](https://tadpole.grand-challenge.org/) challenge. You could find the spreadsheets on `Study Data -> Test Data -> TADPOLE Challenge Data` on [IDA](https://ida.loni.usc.edu/) website.

For AIBL and MACC data, we ran FreeSurfer v6.0 `recon-all` to process T1 MRI data.

----

## Replication

### 1. Setup

Please setup the conda environment before running replication code.

```
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE

conda env create -f replication/config/CBIG_gcVAE_python_env.yml
```

In CBIG, we run experiments via PBS on HPC; if you want to run the code on local workstations, you need to run *every shell script* under `replciations/scripts` folder instead of running the wrapper script `CBIG_gcVAE_replication.sh`; if you want to run the code via other job schedulers, you may need to modify the *submission code* accordingly.


### 2. Run replication experiments

We provied code for replicating three main analyses in [An et al., 2022](https://doi.org/10.1016/j.neuroimage.2022.119570). `unmatch2match` refers to experiments results in section 3.1, 3.2 and 3.3; `sample_size` refers to the experiments results in section 3.4.1;`match2unmatch` refers to experiments results in section 3.4.4. These prefixs (`unmatch2match`, `sample_size` and `match2unmatch`) are used in replication scripts.


The wrapper script `replciations/CBIG_gcVAE_replication.sh` contains all steps to replicate the results in [An et al., 2022](https://doi.org/10.1016/j.neuroimage.2022.119570). For CBIG users, you just need to run the following commands.

```
ssh headnode

cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE

bash replication/CBIG_gcVAE_replication.sh
```

All the results could be found in `results` folder.


### 3. Compare with reference results

Once you have finished the `CBIG_gcVAE_replication.sh`, you could run the following command to check whether you have successfully replicated results.

```
bash replication/scripts/CBIG_gcVAE_replica_check_results.sh
```


### 4. Clean up

You may want to run the following commands to clean up all generated files.

```
rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE/data/match2unmatch

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE/data/unmatch2match

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE/data/sample_size

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE/data/splits

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE/data/raw_data

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE/checkpoints

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE/results

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE/job_logs
```
----

## Bugs and Questions
Please contact Lijun An at anlijuncn@gmail.com, Pansheng Chen at chenpansheng@gmail.com and Chen Zhang at chenzhangsutd@gmail.com
