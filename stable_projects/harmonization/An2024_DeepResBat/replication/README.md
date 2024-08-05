# Replication of DeepResBat harmonization

## References

-   An, L., Zhang, C., Wulan, N., Zhang, S., Chen, P., Ji, F., Ng, KK., Chen, C.,Zhou, J., Yeo, B.T., 2024. [DeepResBat: deep residual batch harmonization accounting for covariate distribution differences](https://doi.org/10.1101/2024.01.18.574145), BioRxiv

---

## Data

If you do not have access to [ADNI](http://adni.loni.usc.edu/), [AIBL](https://aibl.csiro.au/) and [MACC](http://www.macc.sg/) data, you need to apply for these datasets to replicate the results in [An et al., 2024](https://doi.org/10.1101/2024.01.18.574145).

For ADNI data, we used the spreadsheets by [TADPOLE](https://tadpole.grand-challenge.org/) challenge. You could find the spreadsheets on `Study Data -> Test Data -> TADPOLE Challenge Data` on [IDA](https://ida.loni.usc.edu/) website.

For AIBL and MACC data, we ran FreeSurfer v6.0 `recon-all` to process T1 MRI data.

For data matching, please refer to our [released code](https://github.com/ThomasYeoLab/Standalone_An2022_gcVAE/tree/master/matching).

---

## Replication

### 1. Setup

Our experiments were run on GPU servers with `CentOS 7.9`, and the Linux kernel version is `3.10.0`. Please setup the conda environment (e.g., Python version) before running replication code.

```
cd $CBIG_CODE_DIR/stable_projects/harmonization/An2024_DeepResBat

conda env create -f replication/config/CBIG_DeepResBat_python_env.yml
```

In CBIG, we run experiments via PBS on HPC; if you want to run the code on local workstations, you need to run _every shell script_ under `replciations/scripts` folder instead of running the wrapper script `CBIG_DeepResBat_replication.sh`; if you want to run the code via other job schedulers, you may need to modify the _submission code_ accordingly.

### 2. Run replication experiments

We provied code for replicating main analyses (`dataset prediction` and `association analysis`) in [An et al., 2024](https://doi.org/10.1101/2024.01.18.574145).

The wrapper script `replciations/CBIG_DeepResBat_replication.sh` contains all steps to replicate the results in [An et al., 2024](https://doi.org/10.1101/2024.01.18.574145). For CBIG users, you just need to run the following commands.

```
ssh headnode

cd $CBIG_CODE_DIR/stable_projects/harmonization/An2024_DeepResBat

bash replication/CBIG_DeepResBat_replication.sh
```

All the results could be found in `results` folder.

### 3. Compare with reference results

Once you have finished the `CBIG_DeepResBat_replication.sh`, you could run the following command to check whether you have successfully replicated results.

```
python -m utils.check_reference_results
```

### 4. Clean up

You may want to run the following script to clean up all generated files.

```
bash replication/scripts/CBIG_DeepResBat_replica_clean.sh
```

---

## Bugs and Questions

Please contact Lijun An at anlijuncn@gmail.com, and Chen Zhang at chenzhangsutd@gmail.com
