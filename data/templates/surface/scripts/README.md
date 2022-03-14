The `very_inflated` surface still have some foldings and cannot show the insula very well. We inflated this surface a bit more to get a `super_inflated` surface. 

Note that this `super_inflated` surface is not standard and for debug only. 

### Code

`CBIG_super_inflated.sh` generate a `super_inflated` surface from `very_inflated` surface. This script uses `wb_command -surface-inflation`. 

To run the code:

```
mkdir <output_dir>
sh CBIG_super_inflated.sh <output_dir>
```

See [https://www.humanconnectome.org/software/workbench-command/-surface-inflation](https://www.humanconnectome.org/software/workbench-command/-surface-inflation) for more details.

Output surface files are already included in repo: 

```
$CBIG_CODE_DIR/data/templates/surface/fs_LR_32k/fsaverage.L.super_inflated.32k_fs_LR.surf.gii
$CBIG_CODE_DIR/data/templates/surface/fs_LR_32k/fsaverage.R.super_inflated.32k_fs_LR.surf.gii
```

To convert the `.surf.gii` files to freesurfer format, run the following code in matlab:

```
lh = gifti(<your_lh.surf.gii>);
rh = gifti(<your_rh.surf.gii>);
[~, lh_faces] = read_surf(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', 'fs_LR_32k', 'surf', 'lh.very_inflated'));
[~, rh_faces] = read_surf(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', 'fs_LR_32k', 'surf', 'rh.very_inflated'));
write_surf(lh.vertices, lh_faces, <your_lh_output>);
write_surf(rh.vertices, rh_faces, <your_lh_output>);
```

Freesurfer format files are also included:

```
$CBIG_CODE_DIR/data/templates/surface/fs_LR_32k/surf/lh.super_inflated
$CBIG_CODE_DIR/data/templates/surface/fs_LR_32k/surf/rh.super_inflated
```

### Usage

#### Workbench

Load `fsaverage.L.super_inflated.32k_fs_LR.surf.gii` and `fsaverage.R.super_inflated.32k_fs_LR.surf.gii` in wb_view.

#### Matlab

```
CBIG_DrawSurfaceMaps_fslr(lh_labels, rh_labels, 'fs_LR_32k', 'super_inflated')
```


