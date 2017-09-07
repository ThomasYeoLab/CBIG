# Visualize data as annotation on fs_LR surface space

The functions in this folder can visualize data on fs_LR surface space with a better colorscale than what the command line arguments from FreeView allow.

Here is the procedure:
- The surface data on both hemispheres are discretized into discrete bins of values.
- The discretized data are mapped to an annotation map with each label corresponds to a bin of value. The discretized bins of values correpond to discrete values in a colorscale. Thus, the colortable of the annotation map corresponds to the colorscale.
- The resulting annotation maps are automatically opened in FreeView and have their views captured.


# Usage

* **Preparation**

A sample data in MNI space is provided in `./sample` folder. The data can be loaded and projected to fsaverage surface as follows:

```
% Load sample data in volume
x = MRIread('sample/sample_vol.nii.gz');

% Project data in the volume to fsaverage surface
[lh_data, rh_data] = CBIG_ProjectMNI2fsaverage(x, 'fsaverage');
```

* **Visualize data on fs_LR surface**

```
CBIG_DrawFsaverageDataAsFsLRAnnot(lh_data, rh_data, '/data/users/ngohgia/tmp/visualization', 'component', 'parula', 28, 1e-5, 5e-5);
```
Various views of the brains and compiled images of the views are saved under `/data/users/ngohgia/tmp/visualization`, here is one example image:

![visualization_in_fslr](sample/sample_image.png)

* **Visualize data on fs_LR surface with aparc_annot parcellation outlines**

```
CBIG_DrawFsaverageDataAsFsLRAnnot_withParcelOutlines(lh_data, rh_data, '/data/users/ngohgia/tmp/visualization', 'component', 'parula', 'aparc_annot', 28, 1e-5, 5e-5)
```
Here is one example image:

![visualization_in_fslr](sample/sample_image_with_parcellation.png)

