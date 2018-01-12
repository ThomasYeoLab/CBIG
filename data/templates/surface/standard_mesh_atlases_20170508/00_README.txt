The following are the fs_LR standard spheres:

fsaverage.L_LR.spherical_std.164k_fs_LR.surf.gii
L.sphere.59k_fs_LR.surf.gii
L.sphere.32k_fs_LR.surf.gii
fsaverage.R_LR.spherical_std.164k_fs_LR.surf.gii
R.sphere.59k_fs_LR.surf.gii
R.sphere.32k_fs_LR.surf.gii

Importantly, note that "fs_LR" spheres have correspondence between the vertices of the L and R hemispheres.  Some names start with "fsaverage" for historical reasons (and are kept for compatibility purposes; see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3432236); they are not the same as the fsaverage spheres discussed in the next section.

The resample_fsaverage directory contains files to allow resampling between fs_LR and fsaverage with minimal preparatory steps.  See FAQ 9 here:

https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ#HCPUsersFAQ-9.HowdoImapdatabetweenFreeSurferandHCP?

The resample_fsaverage directory contains several resolutions of fsaverage atlas spheres (from FreeSurfer, but converted to GIFTI format)
fsaverage_std_sphere.{L,R}.164k_fsavg_{L,R}.surf.gii (identical to the spheres in (1) above)
fsaverage6_std_sphere.{L,R}.41k_fsavg_{L,R}.surf.gii
fsaverage5_std_sphere.{L,R}.10k_fsavg_{L,R}.surf.gii
fsaverage4_std_sphere.{L,R}.3k_fsavg_{L,R}.surf.gii

and fs_LR spheres deformed to register to the fsaverage atlas
fs_LR-deformed_to-fsaverage.{L,R}.sphere.{164,59,32}k_fs_LR.surf.gii

In order to assist with resampling group average data, the resample_fsaverage directory also contains data files consisting of the group average of the vertex areas ("va_avg") from midthickness surfaces from many HCP subjects (i.e., the various *midthickness_va_avg*.shape.gii files).  These are used with the -area-metrics option of the wb_command -*-resample commands, when appropriate.  (See FAQ mentioned above for some examples).

Tim Coalson, Matt Glasser, Mike Harms, David Van Essen (22 September, 2016)

