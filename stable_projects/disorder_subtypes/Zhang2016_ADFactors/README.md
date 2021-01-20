## Reference

Xiuming Zhang, Elizabeth C. Mormino, Nanbo Sun, Reisa A. Sperling, Mert R. Sabuncu, B. T. Thomas Yeo. [**Bayesian model reveals latent atrophy factors with dissociable cognitive trajectories in Alzheimer's disease**](http://dx.doi.org/10.1073/pnas.1611073113). Proceedings of the National Academy of Sciences of the USA, 113(42):E6535-44, 2016.

----

## Background

Alzheimer’s disease (AD) affects 10% of the elderly population. The disease remains poorly understood with no cure. The main symptom is memory loss, but other symptoms might include impaired executive function (ability to plan and accomplish goals; e.g., grocery shopping). The severity of behavioral symptoms and brain atrophy (gray matter loss) can vary widely across patients. This variability complicates diagnosis, treatment, and prevention. A mathematical model reveals distinct brain atrophy patterns, explaining variation in gray matter loss among AD dementia patients. The atrophy patterns can also explain variation in memory and executive function decline among dementia patients and at-risk nondemented participants. This model can potentially be applied to understand brain disorders with varying symptoms, including autism and schizophrenia.

----

## Data Release 

Factor compositions of the ADNI subjects involved in this study are released as a spreadsheet in the `ADNIDataRelease` folder. Specifically, the spreadsheet includes 810 ADNI1-enrolled subjects’ baseline factor compositions and their factor compositions 24 months after baseline (N = 560).

Folder `ADNIDataRelease/inferNew/files_LDA/model_K3` contain the probabilistic atrophy maps of the atrophy factors for the three-factor model. Maps for the two- and four-factor models are also released (folders `model_K2` and `model_K4`). You can also view or download all of these atrophy maps online at [http://neurovault.org/collections/1917/](http://neurovault.org/collections/1917/).

In folder `ADNIDataRelease/inferNew/`, we provide instructions and code for inferring factor compositions of your own subjects.

----

## Code Release

The remaining folders offer our full code, allowing you to start from scratch on possibly other brain disorders.

Note that this project uses generic functions from other folders, which may be updated over time. To download the exact code utilized for the PNAS paper, you can either

- visit this link:
[https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.18.1-Update_stable_project_unit_test](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.18.1-Update_stable_project_unit_test)

or

- run the following command, if you have Git installed
 
```
git checkout -b Zhang2016_ADFactors v0.18.1-Update_stable_project_unit_test
```

----

## Updates

- Release v0.18.1 (20/01/2021): update unit test to accommodate to the new HPC.

- Release v0.15.3 (16/10/2019): update reference

- Release v0.1.3 (26/04/2017): fixing a bug in wrapper CBIG_LDA_est.sh

- Release v0.1.2 (03/04/2017): removing hardcoded image dimensions

- Release v0.1.1 (11/06/2016): adding README for each internal folder

----

## Bugs and Questions

Please contact xiuming Zhang at xiuming6zhang@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
