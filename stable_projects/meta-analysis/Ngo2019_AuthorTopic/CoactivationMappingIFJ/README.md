# Co-activation Mapping of the left Inferior Frontal Junction (IFJ)

----

## Reference
Gia H. Ngo, Simon B. Eickhoff, Minh Nguyen, Peter T. Fox,  R. Nathan Spreng, B. T. Thomas Yeo. [Beyond Consensus: Embracing Heterogeneity in Curated Neuroimaging Meta-Analysis](https://www.biorxiv.org/content/early/2017/06/19/149567). BioRxiv preprint

Gia H. Ngo, Simon B. Eickhoff, Peter T. Fox, B. T. Thomas Yeo. [Collapsed variational Bayesian inference of the author-topic model: application to large-scale coordinate-based meta-analysis](https://ieeexplore.ieee.org/document/7552332). PRNI2016.

----

# Input text file format
```
//Chikazoe, 2009: Stop vs Go
//Task: 1
-8  28 -11
-10  79 -23

//Zago, 2001: Compute vs Read
//Task: 2
22  24  34
-12   8  51
22 -32 12

...

//Wagner, 2001: 4 Target > 2 Target
//Task: 323
12  13  23
-23   5  35
42 -20 -12
```

- Each experiment has the following format:
  - First line: starts with "//", followed by the identifier of the experiment
  - Second line: starts with "//Task: ", a **UNIQUE** identifier or each experiment, such as the index of the experiment (1, 2, 3, ...) in the dataset. When applying the author-topic model to co-activation mapping, we consider each experiment to employ its own unique task (see the input format when applying the author topic model to estimate functional sub-domains [here](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/meta-analysis/Ngo2019_AuthorTopic/SelfGeneratedThought/MNI152_ActivationCoordinates))
  - Third line onward: coordinates of activation foci reported in MNI152 space.
- Information of one experiment is separated from another by a single empty line

Please note that the activation coordinates used in our IFJ co-activation mapping experiment of the IFJ can only be made available after a data sharing agreement has been obtained with [BrainMap](http://brainmap.org).
We are happy to share the IFJ dataset with any user with the agreement.

---

## Setup
Create a new workspace in Bash (please change this to your preferred directory):
```
mkdir /Work/IFJ_CoactivationMapping
cd /Work/IFJ_CoactivationMapping
matlab
```

---

## Add paths to functions

```
% add paths to functions specific to author-topic model
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'preprocessing'));
addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'inference'));
addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'BIC'));
addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'visualization'));
```

---

## Pre-process MNI152 coordinates of activation foci in raw text input data

Assuming the input data for co-activation mapping of the IFJ is saved at `/Work/data/IFJ_Coordinates.txt` in the text format described at the top.

From Matlab, run the following commands to convert data in the text file to suitable Matlab format for the Collapsed Variational Bayes (CVB) algorithm:
```
textFilePath = '/Work/data/IFJ_Coordinates.txt';
dataDirPath = fullfile(pwd, 'data');
dataFileName = 'IFJ_CVBData.mat';
system(['mkdir -p ' dataDirPath]);

CBIG_AuthorTopic_GenerateCVBDataFromText(textFilePath, dataDirPath, dataFileName);
```

The command above would do the following:
- Save data of the experiments (MNI152 coordinates and task label) at `<dataDirPath>/ExperimentsData.mat`.
- Write the activation foci into brain images under the folder `<dataDirPath>/ActivationVolumes`.
- Perform binary smoothing of the brain activation images and save the resulting images under `<dataDirPath>/BinarySmoothedVolumes`.
- Combine all smoothed brain images into a brain mask of activation across all experiments and save at `<dataDirPath>/mask/expMask.nii.gz`. This brain mask is necessary for computing Bayesian Information Criterion (BIC) measure of the model estimates.
- Save input data for the CVB algorithm at `<dataDirPath>/dataFileName`, i.e. `<dataDirPath>/IFJ_CVBData.mat`.

---

## Inference
To estimate the model parameters with K = 1 to 5 co-activation patterns using 1000 random re-initialization for each K, run the following commands:
```
alpha = 100;
eta = 0.01;
doSmoothing = 1;
workDir = pwd;
cvbData = fullfile(dataDirPath, dataFileName);
seeds = 1:1000;

for K = 1:5
  for seed = seeds
    CBIG_AuthorTopic_RunInference(seed, K, alpha, eta, doSmoothing, workDir, cvbData);
  end
end
```
The 1000 estimates from the 1000 re-initializations for each value of K would be stored under `<workDir>/outputs/K<K>`

---


## Get the best model estimates
To get the estimate with the highest variational log likelihood from the 1000 estimates at each K = 1 to 5,  run the following commands in Matlab:
```
outputsDir = fullfile(workDir, 'outputs');

for K = 1:5
  CBIG_AuthorTopic_ComputeBestSolution(outputsDir, K, seeds, alpha, eta);
end

```
The best solutions would be saved at `<workDir>/outputs/bestSolution/BestSolution_K001.mat` to `<workDir>/outputs/bestSolution/BestSolution_K005.mat`

---


## Visualize the co-activation patterns
To visualize the spatial map of each co-activation pattern (Pr(voxel | co-activation pattern)) on the brain surface for K = 3 components, run the following commands:
```
figuresDir = fullfile(workDir, 'figures');
inputFile = fullfile(workDir, 'outputs', 'bestSolution', 'BestSolution_K003.mat');

CBIG_AuthorTopic_VisualizeComponentsOnBrainSurface(inputFile, figuresDir);
```
This would produce the surface visualization of the K=3 co-activation patterns from the best estimate produced from the previous steps. The images are saved under `figuresDir`.

---

## Compute Bayesian Information Criterion (BIC)

To compute the BIC measures of the best model parameter estimates for K = 1 to 5, run the following commands:

```
maskPath = fullfile(dataDirPath, 'mask', 'expMask.nii.gz');
bestSolutionDir = fullfile(workDir, 'outputs', 'bestSolution');
bicDir = fullfile(workDir, 'BIC');

CBIG_AuthorTopic_ComputeBIC(1:5, maskPath, bestSolutionDir, bicDir);
```
