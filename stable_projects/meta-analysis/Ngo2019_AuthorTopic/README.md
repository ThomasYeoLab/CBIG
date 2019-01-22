## Reference
Gia H. Ngo, Simon B. Eickhoff, Minh Nguyen, Peter T. Fox,  R. Nathan Spreng, B. T. Thomas Yeo. [Beyond Consensus: Embracing Heterogeneity in Neuroimaging Meta-Analysis](https://www.biorxiv.org/content/early/2017/06/19/149567). BioRxiv preprint

----

## Data Release
The activation coordinates used in our experiment for discovering cognitive components of self-generated thought are included under `SelfGeneratedThought/MNI152_ActivationCoordinates`. Explanation of the data format is included the folder's README.

----

## Overview

This project provides all functions necessary to perform a coordinate-based meta-analysis using the author-topic model, such as to estimate cognitive components of a functional domain like self-generated thought (Figure 1) or co-activation patterns of a seed region like the inferior frontal junction (Figure 2).

1. `utilities` folder contains function to perform inference using Collapsed Variational Bayes (CVB) algorithm, compute Bayesian Information Criterion of the estimate, and visualize the estimates.
2. `SelfGeneratedThought` folder contains data and instruction necessary to reproduce the cognitive components of self-generated thought estimated in our paper.
3. `config` folder contains example configuration files compatbile with this project.
4. `unit_tests` folder contains a unit test for this project.

<img src="readme_figures/self_generated_thought.png" width="100%"/>

<img src="readme_figures/ifj.png" width="100%" style="padding-top: 20px"/>

----

## Co-activation Mapping
The author-topic model can be used for co-activation mapping (such as to discover co-activation patterns of the IFJ in our paper (Figure 2 above)).
Users can apply the same data preparation and procedure found under `SelfGeneratedThought/README.md`.
Please note that the acitvatation coordinates used in our IFJ co-activation mapping experiemnt can only be made available after a data sharing agreement has been obtained with [BrainMap](http://brainmap.org).
We are happy to share the IFJ dataset with any user with the agreement.

----

## Configuration
At the time of release, this project was implemented with specific versions of other third-party softwares for computation and visualization, in particular:

* [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) 5.30,
* [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation) 5.0.8,
* [AFNI](https://afni.nimh.nih.gov/pub/dist/tgz/AFNI_ARCHIVE) 16.0.0
* Matlab R2014a

Please see [CBIG repository setup instruction](https://github.com/ThomasYeoLab/CBIG/tree/master/setup) to make your local environment compatible with CBIG repository. The configuration file `setup/CBIG_Ngo2019AuthorTopic_config.sh` or `setup/CBIG_Ngo2019AuthorTopic_config.csh` are example configuration files compatible with this project.

----

## Code Release
### Download stand-alone repository

Since the whole CBIG repository is too big, we provide a stand-alone version of only this project and its dependencies. To download the stand-alone repository, visit this link:
[https://github.com/ThomasYeoLab/Standalone_Ngo2019_AuthorTopic](https://github.com/ThomasYeoLab/Standalone_Ngo2019_AuthorTopic)


### Download whole repository

Except for this project, if you want to use the code for other stable projects from our lab as well, you need to download the whole repository.

- To download the version of the code that was last tested, you can either

  - visit this link:
  [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.8.0-Ngo2019_AuthorTopic](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.8.0-Ngo2019_AuthorTopic)

  or

  - run the following command, if you have Git installed

  ```
  git checkout -b Ngo2019_AuthorTopic v0.8.0-Ngo2019_AuthorTopic
  ```

----

## Bugs and Questions

Please contact Gia at ngohoanggia@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
