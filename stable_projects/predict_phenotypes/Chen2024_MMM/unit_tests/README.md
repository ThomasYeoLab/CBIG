# Unit test for Multilayer Meta-matching (MMM)

The unit test of MM contains two parts:
1. Kernel ridge regression (MATLAB)
2. Deep neural network (Python)

For MATLAB code, we have tested in R2018b.

----

## References
+ Chen, P., An, L., Wulan, N., Zhang, C., Zhang, S., Ooi, L. Q. R., ... & Yeo, B. T. (2023). [**Multilayer meta-matching: translating phenotypic prediction models from multiple datasets to small data**](https://www.biorxiv.org/content/10.1101/2023.12.05.569848v1.abstract). bioRxiv, 2023-12.

----

## Data
The data are synthetic data.

----

## Set up
Please setup the conda environment before running unit tests.
```bash
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM
conda env create -f CBIG_MMM_python_env.yml
conda activate CBIG_Chen2024
```

----
## Run

### Classical KRR
Run the following code in matlab:
```MATLAB
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
base_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'Chen2024_MMM');
unittest_dir = fullfile(base_dir, 'unit_tests');
cd(unittest_dir)

runtests('CBIG_MMM_unit_test')
```

### Transfer learning & Meta-matching methods (Python)
Run the following code in terminal (need gpu to run):
```sh
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/unit_tests/
sh ./CBIG_MMM_unit_test_python.sh
```

----

## Bugs and Questions
Please contact Pansheng Chen at chenpansheng@gmail.com, Lijun An at anlijun.cn@gmail.com, and Chen Zhang at chenzhangsutd@gmail.com.