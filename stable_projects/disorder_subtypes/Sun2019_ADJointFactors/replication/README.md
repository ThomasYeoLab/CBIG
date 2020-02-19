This folder contains scripts about how to replicate results of this project. Notice that all filenames and directories in these scripts **only work for CBIG lab**.

----

## Data
The scripts uses the Alzheimer's Disease Neuroimaging Initiative (ADNI) data `${CBIG_MMLDA_ANDI_DOC_DIR}`. The SPM VBM results of the T1 images are big (2T), you can find it here
```
${CBIG_MMLDA_ADNI_DIR}/Sun2019_SPMVBM
```

----

## Code

There are four main scripts for replication:

1. Get RID, GENDER, AGE, DX, ICV, gmVol and behavioral scores of subjects
To run the following code you need to specify the output directory `subinfo_dir`

```matlab
CBIG_MMLDA_get_subinfo_wrapper(subinfo_dir)
```
To compare results, you can compare with folder

```
${CBIG_REPDATA_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA/data
```

2. Convert grey matter density maps (from VBM) and behavior scores (from above) to documents
To run the following code you need to specify the `subinfo_dir` from step 1 and output directory `doc_dir`

```matlab
CBIG_MMLDA_brain_behavior_to_doc_wrapper(subinfo_dir, doc_dir)
```

To compare results, you can compare with folder

```
${CBIG_REPDATA_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA/results/BrainBehavior2doc
```

3. Run MMLDA estimation on ADNI2 baseline AD and ADNI1 baseline AD subjects
To run the following code you need to specify the `doc_dir` from step2 and output directory `MMLDA_dir`

```bash
./CBIG_MMLDA_runMMLDA_est_wrapper.sh $doc_dir $MMLDA_dir
```

To compare results, you can compare with folder

```
${CBIG_REPDATA_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA/results/estimation
```

4. Visualize factors
To run the following code, you have to make sure that the resolution of VNC is 2000x1000 and "1 big 3 small" mode. Otherwise, taking screenshots of slices of MNI volumn will not work properly.  
To run the following code you need to specify the `mmlda_dir` from step3 and output directory `visualize_dir`

```matlab
CBIG_MMLDA_visualize_factors_wrapper(mmlda_dir, visualize_dir)
```

To compare results, you can compare with folder

```
${CBIG_REPDATA_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA/results/visualizeFactors/
```

5. Run MMLDA inference on ADNI2 baseline MCI, CN, ADNI2 m12 ALL, ADNI2_PETtau
To run the following code you need to specify the `doc_dir` from step2, `visualize_dir` from step4, and output directory `out_dir`

```bash
./CBIG_MMLDA_runMMLDA_inf_wrapper.sh $doc_dir $visualize_dir $out_dir
```

To compare results, you can compare with folder

```
${CBIG_REPDATA_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA/results/inference
```

----
