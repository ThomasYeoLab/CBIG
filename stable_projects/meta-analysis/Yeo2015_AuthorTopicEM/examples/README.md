## Overview

These are examples of performing a small-size meta-analysis using the author-topic modelling of coordinated-based brain activation data.
See (Yeo et al., Cerebral Cortex, 2015) and (Bertolero, Yeo, D'Esposito, PNAS, 2015) for details.

Two examples of parameters estimate using Expectation-Maximization (EM) inference algorithm are provided:
1. `CBIG_AuthorTopicEM_RunEMwithRandInitExample.m` performs EM inference using random initialization
2. `CBIG_AuthorTopicEM_RunEMwithGibbsInitExample.m` performs EM inference using initialization from Gibbs sampling.

The examples were performed on a dataset of 7 tasks engaging self-generated thought, generously provided from the following studies:
- Spreng, R.N., Mar, R.A. and Kim, A.S., 2009. The common neural basis of autobiographical memory, prospection, navigation, theory of mind, and the default mode: a quantitative meta-analysis. Journal of Cognitive Neuroscience.
- Mar, R.A., 2011. The neural bases of social cognition and story comprehension. Annual Review of Psychology.
- Sevinc, G. and Spreng, R.N., 2014. Contextual and perceptual brain processes underlying moral cognition: a quantitative meta-analysis of moral reasoning and moral emotions. PloS One.

In Yeo et al., 2015, we conducted a meta-analysis on the [BrainMap](http://brainmap.org) dataset of 17,000+ neuroimaging experiments, which requires a data sharing agreement to obtain.
Model parameters were estimated using the EM algorithm with initialization from Gibbs sampling.

Please see Ngo et al. 2019 for a new inference algorithm - Collapsed Variational Bayes inference - that is more robust to model's hyperparameters and thus more suitable for common small to common-sized datasets (few hundreds of experiments) and does not not Gibbs sampling for initialization.

## Reference
- BT Thomas Yeo, et al. "Functional specialization and flexibility in human association cortex." Cerebral Cortex, 2015.
- Bertolero, Maxwell A., BT Thomas Yeo, and Mark D'Esposito. "The modular and integrative functional architecture of the human brain." Proceedings of the National Academy of Sciences, 2015.
- Gia H. Ngo, Simon B. Eickhoff, Minh Nguyen, Gunes Sevinc, Peter T. Fox,  R. Nathan Spreng, B. T. Thomas Yeo. Beyond Consensus: Embracing Heterogeneity in Curated Neuroimaging Meta-Analysis. Neuroimage, 2019.
