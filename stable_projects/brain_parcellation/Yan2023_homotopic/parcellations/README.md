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

Note that the ROIs are identical across these three different versions, but only differ in their assignment to different networks.


Homotopic Correspondence
=====================
As mentioned in the paper, the parcellation boundary for the left and right hemispheres do not overlap perfectly. However, each homotopic parcel pair typically occupy roughly symmetric locations across the hemispheres.

Let's take the 400-region parcellation as an example. On the left hemisphere, the parcels are ordered from 1 ~ 200; on the right hemispheres, the parcels are ordered from 201 ~ 400. Parcel 1 correspondes to parcel 201 and so on. Note that for the FreeSurfer version, since there are separate ```.annot``` files for each hemispheres, the parcels are numbered from 1 ~ 200 on within each file.


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
