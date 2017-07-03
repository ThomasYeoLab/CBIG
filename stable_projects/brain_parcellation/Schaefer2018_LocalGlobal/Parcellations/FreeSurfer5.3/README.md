Background
==========
Resting state fMRI data from 1489 subjects were registered using surface-based alignment. A gradient weighted markov random field approach was employed to identify cortical areas.

Information about Downloads
===========================
There are three folders corresponding to the "fsaverage", "fsaverage5" and "fsaverage6" surface space. "fsaverage" contains the high resolution version of the parcellation, while "fsaverage6" and "fsaverage5" contain lower resolution versions of the parcellation. The parcellations were computed in "fsaverage6" space and sampled to "fsaverage5" and "fsaverage". Each parcel was matched to a corresponding network in the 7 and 17 network parcellation by Yeo et al. Please notice that all labels were created in freesurfer version 5.3 and are maybe not be fully functional in other versions of freesurfer.

The structure of each folder follows that of a preprocessed freesurfer subject. In particular, "fsaverage/label/", "fsaverage5/label/", "fsaverage6/label/" contain all the parcellation and confidence files. For example, "fsaverage/label/lh.Schaefer2018_400Parcels_7Networks_order.annot" is the 400 areas parcellation on the left hemisphere ordered and colored according to Yeo et al. 2011.

Example Usage
=============
1) Make sure freesurfer has been set up and configured as usual (http://surfer.nmr.mgh.harvard.edu/fswiki/SetupConfiguration).  
2) In terminal,  
   a) "cd" to unzipped folder containing this README file, fsaverage, fsaverage6, fsaverage5  
   b) freeview -f fsaverage/surf/lh.inflated:annot=fsaverage/label/lh.Schaefer2018_400Parcels_17Networks_order.annot  

References
==========
+ Schaefer A, Kong R, Gordon EM, Zuo XN, Holmes AJ, Eickhoff SB, Yeo BT, (accepted), Local-Global Parcellation of the Human Cerebral Cortex From Intrinsic Functional Connectivity MRI, Cerebral Cortex
+ Yeo BT, Krienen FM, Sepulcre J, Sabuncu MR, Lashkari D, Hollinshead M, Roffman JL, Smoller JW, Zollei L., Polimeni JR, Fischl B, Liu H, Buckner RL (2011) The organization of the human cerebral cortex estimated by functional connectivity. J. Neurophysiol.
