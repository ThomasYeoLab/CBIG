Information about Downloads
===========================
There are three folders corresponding to the `fsaverage`, `fsaverage5` and `fsaverage6` surface space. `fsaverage` contains the highest resolution version of the parcellation, while `fsaverage6` and `fsaverage5` contain lower resolution versions of the parcellation. 

The parcellations were computed in `fsaverage6` space and sampled to `fsaverage5` and `fsaverage`. Please notice that all labels were created in freesurfer version 5.3 and are maybe not be fully functional in other versions of freesurfer.

The structure of each folder follows that of a preprocessed freesurfer subject. In particular, `fsaverage/label`, `fsaverage5/label`, `fsaverage6/label` contain all the parcellation files files. For example, `fsaverage/label/kong17/lh.100Parcels_Kong2022_17Networks.annot` is the 100 areas parcellation on the left hemisphere ordered and colored according to the 17 networks from [Kong et al. 2019](https://pubmed.ncbi.nlm.nih.gov/29878084/).

Example Usage
=============
Make sure freesurfer has been set up and configured as usual (http://surfer.nmr.mgh.harvard.edu/fswiki/SetupConfiguration).  

In terminal,  

1) `cd` to unzipped folder containing this README file, fsaverage, fsaverage6, fsaverage5  

2) run the following command:

```   
freeview -f ./fsaverage/surf/lh.inflated:annot=./fsaverage/label/kong17/lh.100Parcels_Kong2022_17Networks.annot  
```

Bugs and Questions
==================
Please contact Xiaoxuan Yan at xiaoxuan.427@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.