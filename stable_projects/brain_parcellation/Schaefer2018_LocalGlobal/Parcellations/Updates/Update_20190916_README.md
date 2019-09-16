Thanks to feedback from Annie Gupta and her team, we found out that some of the Schaefer parcels have the **wrong label names**, although **the parcel boundaries are correct**. For example, in the 400-region parcellation, all parcels corresponding to Somatomotor Network B are labelled as "auditory" (i.e., `17Networks_LH_SomMotB_Aud_?`), but some of them are obviously not in auditory cortex.

This error arose because of a bug in the algorithm that matches the Schaefer parcellations with the Yeo2011 split components. Note that the Yeo2011 split components came from splitting the original Yeo2011 7/17 networks into spatially contiguous components (https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering). On closer inspection, we also felt that the names of a small number of the Yeo2011 split components were not appropriate.

We have resolved the issues and the updated label names have been released (Release version 0.14.3). See below for more details about label names in Schaefer parcellations.

Parcels' label names
===============
Parcel names are obtained by concatenating the number of networks, hemisphere, network name, component name, and parcel number.
e.g., for `17Networks_LH_DefaultA_pCunPCC_1`:
- `17Networks` means that this parcellation has been matched to the Yeo2011 17 networks
- `LH` means that this parcel is in the left (L) hemisphere (H)
- `DefaultA` is the name of the network that the parcel belongs to
- `pCunPCC` is the component name
- `1` is the component number

Component names are sometimes coarse because they were generated from the Yeo2011 split component parcellation (7/17 networks), where parcels are of relatively low resolution (51/114 regions).

Generating parcels' label names in Schaefer parcellations
==========================================
We developed a new matching algorithm and re-labeled some of the Schaefer parcels. The parcels whose label names have been changed are documented in csv files that provide more details:

+ `${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/Updates/Update_20190916_network_changes.csv`
+ `${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/Updates/Update_20190916_component_name_changes.csv`
+ `${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/Updates/Update_20190916_component_number_changes.csv`

Our updated algorithm does the matching in two steps:
1. The Schaefer parcels are first matched to one of Yeo2011 7/17 networks. Each parcel from Schaefer parcellations is then assigned to a network, e.g., `17Networks_LH_SomMotB`.
2. After that, each Schaefer's parcel is matched to one of the Yeo2011 split components and is assigned a component name, e.g., `17Networks_LH_SomMotB_Cent`.


**IMPORTANT NOTES**
1) When we match the high-resolution Schaefer parcellation (e.g., 1000 parcels) to the Yeo2011 split components, the name might not fully describe each parcel. For example, a Yeo2011 split component might have the name "ParOcc" because it overlaps across parietal (Par) and occipital (Occ) regions. A high-resolution Schaefer parcel might match to the split component and the resulting name is `17Networks_LH_DorsAttnA_ParOcc_?`, but this Schaefer parcel might only be in occipital cortex.

2) Component names are more detailed in the 17 networks than in the 7 networks parcellations.

3) Our lab generally uses the 400-region Schaefer parcellation matched to the Yeo2011 17 networks. 

Details about changes
====

**Changes in parcel ordering**

For most parcels that were incorrectly labeled, only the **component** name has changed. However, for a smaller number of parcels, the **network** name has also been changed. Note that **NOT** all resolutions are affected. For parcellations in which only the component names have changed (but the network names stayed the same), we have kept the ordering of the parcels (in the annot file, colortable, etc) to be the same for consistency sake. However, for parcellations in which the network names have also changed, then we re-ordered the parcels.

**Parcellations whose network names are the same but some component names might have changed**

For example, for the 100-region Schaefer parcellation, parcel `7Networks_LH_DorsAttn_FEF_1` is changed to `7Networks_LH_DorsAttn_PrCv_1`. This is considered as a **change of component name**. 
Parcel `7Networks_LH_DorsAttn_FEF_2` is changed to `7Networks_LH_DorsAttn_FEF_1`. This is considered as a **change of component number**. The table below lists relevant parcellations:

|  Resolution   |  # Parcels whose component name changes  |  # Parcels whose component number changes  |
|  ----  | ----  | ----  |
| 100 Parcels, 7 Networks | 16 | 3 |
| 200 Parcels, 7 Networks  | 29 | 7 |
| 300 Parcels, 7 Networks  | 42 | 6 |
| 400 Parcels, 7 Networks  | 58 | 8 |
| 500 Parcels, 7 Networks  | 67 | 9 |
| 600 Parcels, 7 Networks  | 87 | 13 |
| 700 Parcels, 7 Networks  | 97 | 20 |
| 800 Parcels, 7 Networks  | 112 | 3 |
| 1000 Parcels, 7 Networks  | 123 | 2 |
| 100 Parcels, 17 Networks | 18 | 8 |
| 200 Parcels, 17 Networks  | 36 | 21 |
| 300 Parcels, 17 Networks  | 50 | 33 |
| 400 Parcels, 17 Networks  | 70 | 54 |

**Parcellations whose network names have changed for some parcels**

The network name was changed for a small number of parcels in higher-resolution Schaefer parcellations (usually those near the network boundaries). In this case we will reorder the Schaefer parcellation ordering so that the parcels within the same network will be grouped together. For example, for the 900-region 7-network Schaefer parcellation, `7Networks_LH_DorsAttn_FEF_1` has been changed to `7Networks_LH_Cont_PFCd_1`, the index of this parcel is changed from 205 to 305. The table below lists relevant parcellations:

|  Resolution   |  # Parcels whose network name changes |  # Parcels whose component name changes |  # Parcels whose component number changes  |
|  ----  | ----  | ----  | ----  |
| 900 Parcels, 7 Networks | 1 | 133 | 24 |
| 500 Parcels, 17 Networks  | 1 | 92 | 63 |
| 600 Parcels, 17 Networks  | 2 | 115 | 93 |
| 700 Parcels, 17 Networks  | 2 | 158 | 94 |
| 800 Parcels, 17 Networks  | 1 | 157 | 107 |
| 900 Parcels, 17 Networks  | 4 | 197 | 140 |
| 1000 Parcels, 17 Networks  | 2 | 219 | 154 |

**NOTE ABOUT THE LIMBIC NETWORK**

The Limbic A/B networks in the 17-Network parcellations were displayed in different colors, but the label names did not distinguish between A/B. These labels have now been changed to provide the Limbic A/B distinction.

Mapping parcellation ordering between versions
====

For some users who have been using the old Schaefer parcellations, we provide a `.mat` file that contains the mapping vector between the old and the updated Schaefer parcellations:

+ `${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/Updates/Update_20190916_index_changes.mat`

This `.mat` file is a 20x3 cell matrix, where 20 corresponds to 20 different resolutions. The first column is the parcellation resolution, the second column is the `#parcel x 1` mapping vector, the third column indicates whether the ordering of the parcellation is changed.

If the user would like to generate the mapping vector on their own, we also provide a script:

+ `${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/Code/CBIG_gwMRF_index_trans_btwn2versions.m`

Bugs and Questions
====

Please contact Ruby Kong at roo.cone@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.









