References 
=====================
+ Yan, X., Kong, R., Xue, A., Yang, Q., Orban, C., An, L., Holmes, A.J., Qian, X., Chen, J., Zuo, X.-N., Zhou, J.H., Fortier, M.V., Tan, A.P., Gluckman, P., Chong, Y.S., Meaney, M., Bzdok, D., Eickhoff, S.B., Yeo, B.T.T., 2022. [**Homotopic local-global parcellation of the human cerebral cortex from resting-state functional connectivity**](https://doi.org/10.1101/2022.10.25.513788)

+ Schaefer A, Kong R, Gordon EM, Laumann TO, Zuo XN, Holmes AJ, Eickhoff SB, Yeo BTT. [**Local-Global parcellation of the human cerebral cortex from intrinsic functional connectivity MRI**](http://people.csail.mit.edu/ythomas/publications/2018LocalGlobal-CerebCor.pdf), *Cerebral Cortex*, 29:3095-3114, 2018

+ Yeo BTT, Krienen FM, Sepulcre J, Sabuncu MR, Lashkari D, Hollinshead M, Roffman JL, Smoller JW, Zöllei L, Polimeni JR, Fischl B, Liu H, Buckner RL. [**The organization of the human cerebral cortex revealed by intrinsic functional connectivity**](http://people.csail.mit.edu/ythomas/publications/2011CorticalOrganization-JNeurophysiol.pdf), *J Neurophysiology*, 106(3):1125-1165, 2011

+ Kong, R., Li, J., Orban, C., Sabuncu, M.R., Liu, H., Schaefer, A., Sun, N., Zuo, X.-N., Holmes, A.J., Eickhoff, S.B., 2019. [**Spatial topography of individual-specific cortical networks predicts human cognition, personality, and emotion**](https://pubmed.ncbi.nlm.nih.gov/29878084/) *Cerebral cortex* 29, 2533–2551


Parcellations Release
=====================
This folder contains parcellations are available at multiple resolutions (100 parcels to 1000 parcels). Specifically, there are three subfolders corresponding to three different spaces ```Freesurfer```, ```MNI``` and ```HCP```. 

The parcellations were first computed in Freesurfer ```fsaverage6``` space and sampled to ```fsaverage5``` and ```fsaverage``` space. The parcellations were also projected to HCP ```fslr32k``` and FSL ```MNI``` space. 

For each resolution, there are three different versions where the ROIs are matched to different large-scale RSFC networks:

* VERSION 1. `<#ROIs>Parcels_Yeo2011_7Networks` under `yeo7` folders.
* VERSION 2. `<#ROIs>Parcels_Yeo2011_17Networks` under `yeo17` folders.
* VERSION 3. `<#ROIs>Parcels_Kong2022_17Networks` under `kong17` folders.

We matched each parcel to a specific network based on maximal spatial overlap, i.e., if a parcel has 80% overlap with network A and 20% with network B, it would be assigned to network A.

Note that the vertices enclosed in each ROI are identical across these three different versions, and the ROIs only differ in their assignment to different networks.

To use the parcellations without the trouble of downloading our entire repository, you can just click on this link: [download Yan2023 Parcellations](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Yan2023_homotopic/parcellations)

Homotopic Correspondence
=====================
As mentioned in the paper, the parcellation boundary for the left and right hemispheres do not overlap perfectly. However, each homotopic parcel pair typically occupy roughly symmetric locations across the hemispheres.

Let's take the 400-region parcellation as an example. On the left hemisphere, the parcels are ordered from 1 ~ 200; on the right hemispheres, the parcels are ordered from 201 ~ 400. Parcel 1 correspondes to parcel 201 and so on. Note that for the FreeSurfer version, since there are separate ```.annot``` files for each hemispheres, the parcels are numbered from 1 ~ 200 on within each file.

Also note that, for each homotopic parcel pair, they might have been assigned to different networks due to the asymmetry of large-scale brain networks. For example, in `./HCP/fsLR32k/kong17/100Parcels_Kong2022_17Networks_info.txt`:

> Line 27~28, parcel 14 `17networks_LH_DefaultA_RSC` corresponds to Line 127-128, parcel 64 `17networks_RH_DefaultC_RSC`. This pair of homotopic parcels belong to different sub default mode networks.

Alternative Parcellation Version
========================
An alternative version of our parcellation is available, where the parcels on each hemisphere are ordered simply by network assignment. Therefore, there is no homotopic correspondence between the left and right hemispheres. You should only consider this version if your analysis does not require homotopic correspondence.

Note that in this alternative version, only the parcels on the right hemispheres are re-ordered, as compared the standard version that we provide in this folder.

Download this alternative version here (part of what we have previously released): [download Yan2023 Parcellations (ordered by network within each hemisphere)](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/ThomasYeoLab/CBIG/tree/c874d05cf760999427d7d2be5390b4fea21dd439/stable_projects/brain_parcellation/Yan2023_homotopic/parcellations)

Parcel Naming Convention
========================

Let's look at a few examples.

The first row in `./MNI/kong17/freeview_lut/100Parcels_Kong2022_17Networks_LUT.txt`:

> 1 17networks_LH_DefaultC_PHC 0 0 130 0

indicates that the current 100-level parcellation was matched to Kong2022 networks, parcel 1 is on the left hemisphere, assigned to the DefaultC network and resides roughly on the PHC cortex. `0 0 130` is just the RGB color code for this parcel, and the last trailing `0` doesn't make any difference.

The last row in `./MNI/yeo7/freeview_lut/200Parcels_Yeo2011_7Networks_LUT.txt`:

> 200 7networks_RH_Vis_Striate_2 120 18 134 0

indicates that the current 200-level parcellation was matched to Yeo7 networks, parcel 200 is on the right hemisphere, assigned to the Visual network and resides roughly on the stiate cortex. `120 18 134` is just the RGB color code for this parcel, and the last trailing `0` doesn't make any difference.

For abbreviations used in the parcel names (anatomic locations), refer to the table in the following section `Parcel Name Abbreviations`.

RAS centroid coordinates
=====================
We also provide RAS centroid coordinates of the Yeo 7/17 and Kong 17 parcellations in MNI 1mm and 2mm space under: `MNI/centroid_coordinates`


Parcel Name Abbreviations
=====================
When naming each parcel, we use a set of abbreviations to indicate various anatomical locations. The corresponding full parcel names can be found here:

| Abbreviation | Full parcel Name |
| ---- | ---- |
| Cingm | mid-cingulate |
| ExStr | extrastriate cortex |
| ExStrInf | extra-striate inferior |
| ExStrSup | extra-striate superior |
| FPole | frontal pole |
| FrMed | frontal medial |
| FrOper | frontal operculum |
| Ins | insula |
| IPL | inferior parietal lobule |
| IPS | intraparietal sulcus |
| OFC | orbital frontal cortex |
| ParMed | parietal medial |
| PCC | posterior cingulate cortex |
| pCun | precuneus |
| pCunPCC | precuneus posterior cingulate cortex |
| PFCd | dorsal prefrontal cortex |
| PFCl | lateral prefrontal cortex |
| PFCm | medial prefrontal cortex |
| PFCv | ventral prefrontal cortex |
| PHC | parahippocampal cortex |
| PostC | post central |
| PrC | precentral |
| PrCd | precentral dorsal |
| PrCv | precentral ventral |
| RSC | retrosplenial cortex |
| SPL | superior parietal lobule |
| ST | superior temporal |
| Striate | striate cortex |
| Temp | temporal |
| TempOcc | temporal occipital |
| TempPole | temporal pole |
