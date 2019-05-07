## Probabilistic atrophy maps and cognitive deficits
In our paper, we estimated probabilistic atrophy maps and cognitive deficits of factors on both ADNIGO/2 and ADNI1 independently. Even though they are quite similar, language is assoicated with left temporal atrophy in ADNIGO/2 and assoicated with bilateral temporal atrophy in ADNI1. Therefore, we provided both ADNIGO/2 and ADNI1 results separately.

- ADNIGO/2
    - k=2
        - `MMLDA_files/ADNIGO2/k2/atrophy1.nii.gz`: posterior cortical
        - `MMLDA_files/ADNIGO2/k2/atrophy2.nii.gz`: medial temporal
        - `MMLDA_files/ADNIGO2/k2/behavior1.csv`: executive function
        - `MMLDA_files/ADNIGO2/k2/behavior2.csv`: memory
    - k=3
        - `MMLDA_files/ADNIGO2/k3/atrophy1.nii.gz`: medial temporal
        - `MMLDA_files/ADNIGO2/k3/atrophy2.nii.gz`: posterior cortical
        - `MMLDA_files/ADNIGO2/k3/atrophy3.nii.gz`: lateral temporal
        - `MMLDA_files/ADNIGO2/k3/behavior1.csv`: memory
        - `MMLDA_files/ADNIGO2/k3/behavior2.csv`: executive function
        - `MMLDA_files/ADNIGO2/k3/behavior3.csv`: language
    - k=4
        - `MMLDA_files/ADNIGO2/k4/atrophy1.nii.gz`: medial temporal 2 
        - `MMLDA_files/ADNIGO2/k4/atrophy2.nii.gz`: lateral temporal
        - `MMLDA_files/ADNIGO2/k4/atrophy3.nii.gz`: medial temporal
        - `MMLDA_files/ADNIGO2/k4/atrophy4.nii.gz`: posterior cortical
        - `MMLDA_files/ADNIGO2/k4/behavior1.csv`: mixed
        - `MMLDA_files/ADNIGO2/k4/behavior2.csv`: language
        - `MMLDA_files/ADNIGO2/k4/behavior3.csv`: memory
        - `MMLDA_files/ADNIGO2/k4/behavior4.csv`: executive function

- ADNI1
    - k=3
        - `MMLDA_files/ADNI1/k3/atrophy1.nii.gz`: medial temporal
        - `MMLDA_files/ADNI1/k3/atrophy2.nii.gz`: bilateral temporal
        - `MMLDA_files/ADNI1/k3/atrophy3.nii.gz`: posterior cortical
        - `MMLDA_files/ADNI1/k3/behavior1.nii.gz`: memory
        - `MMLDA_files/ADNI1/k3/behavior2.nii.gz`: language
        - `MMLDA_files/ADNI1/k3/behavior3.nii.gz`: executive function

You can also view or download all of these atrophy maps online at [http://neurovault.org/collections/xxx/](http://neurovault.org/collections/xxx/).

----

## Factor compositions of ADNIGO/2 and ADNI1 participants

Spreadsheet `FactorCompositions_ADNIGO2_895bl-668m12.csv` includes 895 ADNIGO/2-enrolled subjects’ baseline factor compositions and their factor compositions 12 months after baseline (N = 668). Therefore, it has 1563 rows (excluding the header row).  

Spreadsheet `FactorCompositions_ADNI1_777bl.csv` includes 777 ADNI1-enrolled subjects' baseline factor compositions.

### Column descriptions
It has 14 columns and the description of each column: 

1. `RID` (roster ID). This column corresponds to the `RID` column of, e.g., the diagnosis file `DXSUM_PDXCONV_ADNIALL.csv` (downloadable from the ADNI website). It also matches with the last four digits of the `SubjectID` column of `MPRAGEMETA.csv` (downloadable from the ADNI website).
2. `VISCODE2` (visit code). This column corresponds to the `VISCODE2` column of the diagnosis file `DXSUM_PDXCONV_ADNIALL.csv`. Possible values in this spreadsheet are `bl`, standing for baseline (or screening), and `m12`, standing for 12 months after baseline.
3. `ImageUID` (image ID). This column corresponds to the `ImageUID` column of `MPRAGEMETA.csv`. Using `ImageUID`, You can uniquely identify an image within the whole ADNI dataset. 
4. `Phase` (study phase). In this spreadsheet, its value is `ADNIGO` or `ADNI2`.
5. `ScanDate` (scan date). This column corresponds to the `ScanDate` column of `MPRAGEMETA.csv`. Note that this date is different from, but usually very close to, the examination date (the `EXAMDATE` column of the diagnosis file `DXSUM_PDXCONV_ADNIALL.csv`). 
6. `MTL_Memory_Prob`. Medial Temporal Lobe and Memory joint factor probability using both atrophy and behavior.
7. `LateralTemporal_Language_Prob`. Lateral Temporal and Language joint factor probability using both atrophy and behavior.
8. `PosteriorCortical_Executive_Prob`. Posterior Cortical and Executive Function joint factor probability using both atrophy and behavior.
9. `MTL_Prob`. Medial Temporal Lobe factor probability using only atrophy.
10. `LateralTemporal_Prob`. Lateral Temporal factor probability using only atrophy.
11. `PosteriorCortical_Prob`. Posterior Cortical factor probability using only atrophy.
12. `Memory_Prob`. Memory factor probability using only behavior.
13. `Language_Prob`. Language factor probability using only behavior.
14. `Executive_Prob`. Executive function factor probability using only behavior.

----

## SPM VBM ADNIGO/2 and ADNI1 template and mask
We used CAT12 toolbox to create a study specific template 

```
SPM_VBM_files/ADNIGO2/wTemplate_0.nii
SPM_VBM_files/ADNIGO2/wTemplate_1.nii
SPM_VBM_files/ADNIGO2/wTemplate_2.nii
SPM_VBM_files/ADNIGO2/wTemplate_3.nii
SPM_VBM_files/ADNIGO2/wTemplate_4.nii
SPM_VBM_files/ADNIGO2/wTemplate_5.nii
SPM_VBM_files/ADNIGO2/wTemplate_6.nii

SPM_VBM_files/ADNI1/wTemplate_0.nii 
SPM_VBM_files/ADNI1/wTemplate_1.nii 
SPM_VBM_files/ADNI1/wTemplate_2.nii 
SPM_VBM_files/ADNI1/wTemplate_3.nii 
SPM_VBM_files/ADNI1/wTemplate_4.nii 
SPM_VBM_files/ADNI1/wTemplate_5.nii 
SPM_VBM_files/ADNI1/wTemplate_6.nii 
```

and also generate the grey matter mask 

```
SPM_VBM_files/ADNIGO2/ADNI2_bl_gm_mask_MNI2mm.nii.gz
SPM_VBM_files/ADNI1/ADNI1_bl_gm_mask_MNI2mm.nii.gz
```
These files are provided if you want to replicate our results. For more info, see `unit_tests/replicate` folder.

----

## Regression and zscore parameters
After doing the VBM and combining some behavioral scores, we regressed out the 
[age, sex, icv] with respect to cognitive normal (CN) participants. In addition,
we also did zscore with respect to cognitive normal (CN) pariticipants.

In the following mat files, it contains reference (CN) parameters from 
`step2_MMLDA/CBIG_MMLDA_brain_to_zscore_create_refpara.m` and 
`step2_MMLDA/CBIG_MMLDA_behavior_to_zscore_create_refpara.m`.

The variable [ref_mean] [ref_beta] will be used for regression with respect to reference (CN) cohort.
The variable [ref_reg_mean] [ref_reg_std] will be used for z score with respect to reference (CN) cohort. 

```
regression_zscore_paras/ADNI2_bl_brain_refpara.mat
regression_zscore_paras/ADNI2_bl_behavior_refpara.mat
regression_zscore_paras/ADNI1_bl_brain_refpara.mat
regression_zscore_paras/ADNI1_bl_behavior_refpara.mat
```

