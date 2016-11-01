This folder generates Figs. 5, 7, S6, S7 and S8.

This folder utilizes the LME model to explore the relationships between longitudinal cognitive decline rates and atrophy factors.

---

## What Does Each File Do?

1. `CBIG_wrapper.m` shows how to call the functions.
2. `CBIG_fitLME_3t.m` fits a three-factor LME model to the behavioral scores. Ditto for `CBIG_fitLME_2t.m` and `CBIG_fitLME_4t.m`.
3. `CBIG_hypoTest_3t.m` performs statistical tests among the three factors. Ditto for `CBIG_hypoTest_2t.m` and `CBIG_hypoTest_4t.m`.
4. `CBIG_hypoTest_3t_withinFactorAcrossStages.m` performs statistical tests across different disease stages within a given factor.
5. The rest are our p values and raw figures.
