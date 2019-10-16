## References

Wang P, Kong R, Kong XL, Liegeois R, Orban C, Deco G, Van den Heuvel M, Yeo BT. [**Inversion of a Large-Scale Circuit Model Reveals a Cortical Hierarchy in the Dynamic Resting Human Brain**](https://advances.sciencemag.org/content/5/1/eaat7854.full), Science advances, 5(1):eaat7854, 2019

----

## Background

The human brain is a complex multi-scale system. Here, a large-scale model of neural dynamics with region-specific micro-scale properties is considered. A stochastic optimization framework was adopted to invert the model, yielding dramatically better fit to new, out-of-sample functional magnetic resonance imaging data. Without assuming the existence of a hierarchy, the estimated model parameters revealed a large-scale cortical gradient. At one end, sensory-motor cortex possessed strong recurrent connections and excitatory subcortical inputs, consistent with localized processing of external stimuli. At the opposing end, the default network possessed weak recurrent connections and excitatory subcortical inputs, consistent with its role in internal thought. Finally, recurrent connection strength was associated with laminar-specific neuronal cell density (but not cell size) in an independent histological dataset. Overall, this study provides micro-scale insights into a macro-scale cortical hierarchy in the dynamic resting brain. 

----

## Code Release

This folder contains all matlab data, functions and scripts to replicate the result of our project. This code uses dynamic mean-field model to simulate the functional connectivity (FC). Inverse parameter estimation method is applied to fit the model parameter to the empiricial data. 

Our scripts only applicable for FC and SC in Desikan parcellation. Note that Desikan parcellation has 36 ROIs in each hemisphere, but we ignore two ROIs: `unkown` and `corpuscallosum`. In our case, there are 34 ROIs in each hemisphere.

Please read `README.md` in each sub-folder for details on each step.

#### Step 1: Estimation
The sub-folder `step1_estimation/` contains the first step which will fit the model to the empirical data and estimate parameters.

#### Step 2: Simulation
The sub-folder `step2_simulation/` contains the second step which will simulate the FC using the estimated parameters from Step1: Estimation.

#### Step 3: Plot
The sub-folder `step3_plot/` contains all data, functions and scripts to plot the brain map and bar chart of the estimated parameters obtained in Step 1: Estimation.

#### Download
To download the version of the code that is last tested, you can either visit this link: https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.15.3-Update_proj_refs_and_add_KRR_LITE or run the following command, if you have Git installed:
```
git checkout -b Wang2018_MFMem v0.15.3-Update_proj_refs_and_add_KRR_LITE
```

If you don't want to download our entire repository and only want this project `Wang2018_MFMem` only, you can use https://kinolien.github.io/gitzip/ to download this folder and these two extra scripts:
```
$CBIG_CODE_DIR/utilities/matlab/surf/CBIG_ReadNCAvgMesh.m
$CBIG_CODE_DIR/utilities/matlab/figure_utilities/CBIG_DrawSurfaceMapsWithBoundary.m
```

#### Examples
In the sub-folder `examples/`, we provide examples on how to run the code. Moreover, if the users would like to apply our code on their own data, we also provide instructions on how to realize that.

#### Function Library
The sub-folder `lib/` contains all the functions created for this project. 

#### Estimation and Simualtion for the 4 parameter model
The sub-folder `extra_4p` contains all matlab data, functions and scripts to show some results we mentioned in the paper. In this experiment, we use dynamic mean-field model to simulate the functional connectivity (FC) obtained from resting-state fMRI with only 4 parameters. The details of the model can be referred from our paper.

----

## Bugs and Questions

Please contact Dr. Peng Wang (pengwanghome@gmail.com) and Prof. Dr. Thomas Yeo (yeoyeo02@gmail.com).

