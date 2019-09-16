Thanks to feedback from Annie Gupta and her team, we found out that some of the Schaefer parcels have the **wrong label names**, although **the parcel boundaries are correct**. 

This error arose because of a bug in the algorithm that matches the Schaefer parcellations with the Yeo2011 split components. Note that the Yeo2011 split components came from splitting the original Yeo2011 7/17 networks into spatially contiguous components. On closer inspection, we also felt that **the names of a small number of the Yeo2011 split components were not appropriate**.

We have resolved the issues and the updated some label names for Yeo2011 split componnets, and the updated files have been released (Release version 0.14.3). See below for more details about label names.

Label names
===============
Label names are obtained by concatenating the number of networks, hemisphere, network name, component name.
e.g., for `17Networks_LH_DefaultA_pCunPCC`:
- `17Networks` means that it's a component from Yeo2011 7 or 17 networks
- `LH` means that this parcel is in the left (L) hemisphere (H)
- `DefaultA` is the name of the network that the component belongs to
- `pCunPCC` is the component name

Updates
==========================================
We renamed some of the components. The components whose label names have been changed are documented in csv files that provide more details:

+ `${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/Updates/Update_20190916_component_name_changes.csv`
+ `${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/Updates/Update_20190916_network_name_changes.csv`

Bugs and Questions
====

Please contact Ruby Kong at roo.cone@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
