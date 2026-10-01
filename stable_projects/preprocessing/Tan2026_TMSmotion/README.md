# Evaluation of censoring on the quality of parcellations & TMS targets


## References
Tan, T. W. K., Kong, R., Xue, A., Cheng, J., et al. (2026). [Excessive Censoring Degrades Individual-Specific Cortical Parcellations and Personalized TMS Targets](https://doi.org/10.64898/2026.03.09.710457). _bioRxiv_, 2026.2003.2009.710457.

---

## Background
Head motion is a major source of bias in resting-state functional connectivity (FC) estimates. A common strategy is to censor high-motion volumes and discard high-motion runs, but overly stringent censoring may remove meaningful signal along with noise and degrade downstream estimates. In this work, we examine how motion-censoring strategies affect individualized cortical parcellations and personalized transcranial magnetic stimulation (TMS) target selection using precision-fMRI data. We find that higher motion produces greater deviations from individualized low-motion references, but at a given motion level, lenient censoring consistently outperforms strict censoring. The improvement is comparable to doubling scan duration under strict censoring. These results suggest that within the motion range and parcellation/TMS targeting frameworks evaluated here, high-motion sessions may remain useful for individualized parcellation and TMS targeting, potentially reducing the need for re-scanning.

---

## Folder structure
* `select_subID`: Contains Python scripts to select participants that qualify.
* `TMS`: Contains Python scripts to run TMS algorithm.
* `compare`: Contains a Python script to compare similarity to ground-truth.
* `util`: Contains a general utilities script.
* `replication`: Contains the scripts to replicate analyses.
* `examples`: Contains the example data, scripts and reference outputs.
* `unit_tests`: Runs code in `examples` to check code.
* `config`: Contains conda env .yml, and tested config.

---

## Code Release
### Download whole repository
Besides this project, if you want to use the code for other stable projects from our lab as well, you need to download the whole repository.

To download the version of the code that was last tested, you can either: 

- visit this link: [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.41.0-Tan2026_TMSmotion](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.41.0-Tan2026_TMSmotion)

or

- run the following command, if you have Git installed

```
git checkout -b Tan2026_TMSmotion v0.41.0-Tan2026_TMSmotion
```


---

## Usage
### Setup

#### Download Dependency Repository
1. Please download the TMS targeting code using this link: https://github.com/ThomasYeoLab/Tan2026_TMSmotion_target
2. Place downloaded `Tan2026_TMSmotion_target` directory at `$CBIG_CODE_DIR/..` (one level above `$CBIG_CODE_DIR`)

#### Install Python environment
```
cd $CBIG_CODE_DIR/stable_projects/preprocessing/Tan2026_TMSmotion/config
conda env create -f CBIG_TMSmotion_python_env.yml
```

#### Environment variables setup
1. Please run:
```
export TAN2026_TMSMOTION_DIR=$CBIG_CODE_DIR/stable_projects/preprocessing/Tan2026_TMSmotion
export KONG2026_TMS_DIR=$CBIG_CODE_DIR/../Tan2026_TMSmotion_target
```

2. Please `export` correct software versions (using your own paths), detailed in `$CBIG_CODE_DIR/stable_projects/preprocessing/Tan2026_TMSmotion/config/CBIG_TMSmotion_tested_config.sh`

3. To initialize environment variables properly, please restart the terminal.

---

## Updates
* Release v0.41.0 (01/10/2026): Initial release of Tan2026_TMSmotion project
---

## Bugs and questions
Please contact Trevor Tan at [trevortanweikiat97@gmail.com](mailto:trevortanweikiat97@gmail.com) and Thomas Yeo at [yeoyeo02@gmail.com](mailto:yeoyeo02@gmail.com).
