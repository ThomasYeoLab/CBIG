## Reference

Xiuming Zhang, Elizabeth C. Mormino, Nanbo Sun, Reisa A. Sperling, Mert R. Sabuncu, B. T. Thomas Yeo. [**Bayesian model reveals latent atrophy factors with dissociable cognitive trajectories in Alzheimer's disease**](http://dx.doi.org/10.1073/pnas.1611073113). *Proceedings of the National Academy of Sciences of the USA*, 2016.

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
[https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.1-release_Zhang2016_ADFactors](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.1-release_Zhang2016_ADFactors)

or

- run the following command, if you have Git installed
 
```
git checkout -b Zhang2016_ADFactors v0.1-Zhang2016_ADFactors
```

----

## Bugs and Questions

Please contact Xiuming Zhang at (firstname)6(lastname)@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
