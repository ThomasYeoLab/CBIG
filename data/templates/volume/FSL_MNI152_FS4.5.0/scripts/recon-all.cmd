
 mri_convert /autofs/space/lyon_006/pubsw/Linux2-2.3-x86_64/packages/fsl.64bit/4.1.5/data/standard/MNI152_T1_1mm.nii.gz /autofs/cluster/nexus/12/users/ythomas/data/MaiAtlas/FSL_MNI152_FS/mri/orig/001.mgz 

#--------------------------------------------
#@# MotionCor Thu Jul 15 08:40:54 EDT 2010

 cp /autofs/cluster/nexus/12/users/ythomas/data/MaiAtlas/FSL_MNI152_FS/mri/orig/001.mgz /autofs/cluster/nexus/12/users/ythomas/data/MaiAtlas/FSL_MNI152_FS/mri/rawavg.mgz 


 mri_convert /autofs/cluster/nexus/12/users/ythomas/data/MaiAtlas/FSL_MNI152_FS/mri/rawavg.mgz /autofs/cluster/nexus/12/users/ythomas/data/MaiAtlas/FSL_MNI152_FS/mri/orig.mgz --conform 


 mri_add_xform_to_header -c /autofs/cluster/nexus/12/users/ythomas/data/MaiAtlas/FSL_MNI152_FS/mri/transforms/talairach.xfm /autofs/cluster/nexus/12/users/ythomas/data/MaiAtlas/FSL_MNI152_FS/mri/orig.mgz /autofs/cluster/nexus/12/users/ythomas/data/MaiAtlas/FSL_MNI152_FS/mri/orig.mgz 

#--------------------------------------------
#@# Nu Intensity Correction Thu Jul 15 08:41:08 EDT 2010

 mri_nu_correct.mni --i orig.mgz --o nu.mgz --n 2 

#--------------------------------------------
#@# Talairach Thu Jul 15 08:43:05 EDT 2010

 talairach_avi --i nu.mgz --xfm transforms/talairach.auto.xfm 


 cp transforms/talairach.auto.xfm transforms/talairach.xfm 

#--------------------------------------------
#@# Talairach Failure Detection Thu Jul 15 08:43:55 EDT 2010

 talairach_afd -T 0.005 -xfm transforms/talairach.xfm 


 awk -f /usr/local/freesurfer/stable4//bin/extract_talairach_avi_QA.awk /autofs/cluster/nexus/12/users/ythomas/data/MaiAtlas/FSL_MNI152_FS/mri/transforms/talairach_avi.log 

#--------------------------------------------
#@# Intensity Normalization Thu Jul 15 08:43:56 EDT 2010

 mri_normalize -g 1 nu.mgz T1.mgz 

#--------------------------------------------
#@# Skull Stripping Thu Jul 15 08:46:26 EDT 2010

 mri_em_register -skull nu.mgz /usr/local/freesurfer/stable4//average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta 


 mri_watershed -T1 -brain_atlas /usr/local/freesurfer/stable4//average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta T1.mgz brainmask.auto.mgz 


 cp brainmask.auto.mgz brainmask.mgz 

#-------------------------------------
#@# EM Registration Thu Jul 15 09:26:27 EDT 2010

 mri_em_register -mask brainmask.mgz nu.mgz /usr/local/freesurfer/stable4//average/RB_all_2008-03-26.gca transforms/talairach.lta 

#--------------------------------------
#@# CA Normalize Thu Jul 15 09:49:51 EDT 2010

 mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /usr/local/freesurfer/stable4//average/RB_all_2008-03-26.gca transforms/talairach.lta norm.mgz 

#--------------------------------------
#@# CA Reg Thu Jul 15 09:51:28 EDT 2010

 mri_ca_register -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz /usr/local/freesurfer/stable4//average/RB_all_2008-03-26.gca transforms/talairach.m3z 

#--------------------------------------
#@# CA Reg Inv Thu Jul 15 13:38:27 EDT 2010

 mri_ca_register -invert-and-save transforms/talairach.m3z 

#--------------------------------------
#@# Remove Neck Thu Jul 15 13:39:13 EDT 2010

 mri_remove_neck -radius 25 nu.mgz transforms/talairach.m3z /usr/local/freesurfer/stable4//average/RB_all_2008-03-26.gca nu_noneck.mgz 

#--------------------------------------
#@# SkullLTA Thu Jul 15 13:40:19 EDT 2010

 mri_em_register -skull -t transforms/talairach.lta nu_noneck.mgz /usr/local/freesurfer/stable4//average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta 

#--------------------------------------
#@# SubCort Seg Thu Jul 15 14:18:56 EDT 2010

 mri_ca_label -align -nobigventricles norm.mgz transforms/talairach.m3z /usr/local/freesurfer/stable4//average/RB_all_2008-03-26.gca aseg.auto_noCCseg.mgz 


 mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz FSL_MNI152_FS 

#--------------------------------------
#@# Merge ASeg Thu Jul 15 14:43:33 EDT 2010

 cp aseg.auto.mgz aseg.mgz 

#--------------------------------------------
#@# Intensity Normalization2 Thu Jul 15 14:43:33 EDT 2010

 mri_normalize -aseg aseg.mgz -mask brainmask.mgz norm.mgz brain.mgz 

#--------------------------------------------
#@# Mask BFS Thu Jul 15 14:46:46 EDT 2010

 mri_mask -T 5 brain.mgz brainmask.mgz brain.finalsurfs.mgz 

#--------------------------------------------
#@# WM Segmentation Thu Jul 15 14:46:51 EDT 2010

 mri_segment brain.mgz wm.seg.mgz 


 mri_edit_wm_with_aseg -keep-in wm.seg.mgz brain.mgz aseg.mgz wm.asegedit.mgz 


 mri_pretess wm.asegedit.mgz wm norm.mgz wm.mgz 

#--------------------------------------------
#@# Fill Thu Jul 15 14:49:53 EDT 2010

 mri_fill -a ../scripts/ponscc.cut.log -xform transforms/talairach.lta -segmentation aseg.auto_noCCseg.mgz wm.mgz filled.mgz 

#--------------------------------------------
#@# Tessellate lh Thu Jul 15 14:50:52 EDT 2010

 mri_pretess ../mri/filled.mgz 255 ../mri/norm.mgz ../mri/filled-pretess255.mgz 


 mri_tessellate ../mri/filled-pretess255.mgz 255 ../surf/lh.orig.nofix 


 rm -f ../mri/filled-pretess255.mgz 


 mris_extract_main_component ../surf/lh.orig.nofix ../surf/lh.orig.nofix 

#--------------------------------------------
#@# Smooth1 lh Thu Jul 15 14:51:06 EDT 2010

 mris_smooth -nw -seed 1234 ../surf/lh.orig.nofix ../surf/lh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 lh Thu Jul 15 14:51:13 EDT 2010

 mris_inflate -no-save-sulc ../surf/lh.smoothwm.nofix ../surf/lh.inflated.nofix 

#--------------------------------------------
#@# QSphere lh Thu Jul 15 14:52:25 EDT 2010

 mris_sphere -q -seed 1234 ../surf/lh.inflated.nofix ../surf/lh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology lh Thu Jul 15 15:06:02 EDT 2010

 cp ../surf/lh.orig.nofix ../surf/lh.orig 


 cp ../surf/lh.inflated.nofix ../surf/lh.inflated 


 mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 FSL_MNI152_FS lh 


 mris_euler_number ../surf/lh.orig 


 mris_remove_intersection ../surf/lh.orig ../surf/lh.orig 


 rm ../surf/lh.inflated 

#--------------------------------------------
#@# Make Final Surf lh Thu Jul 15 15:53:25 EDT 2010

 mris_make_surfaces -noaparc -mgz -T1 brain.finalsurfs FSL_MNI152_FS lh 

#--------------------------------------------
#@# Surf Volume lh Thu Jul 15 16:55:48 EDT 2010

 mris_calc -o lh.area.mid lh.area add lh.area.pial 


 mris_calc -o lh.area.mid lh.area.mid div 2 


 mris_calc -o lh.volume lh.area.mid mul lh.thickness 

#--------------------------------------------
#@# Smooth2 lh Thu Jul 15 16:55:49 EDT 2010

 mris_smooth -n 3 -nw -seed 1234 ../surf/lh.white ../surf/lh.smoothwm 

#--------------------------------------------
#@# Inflation2 lh Thu Jul 15 16:55:54 EDT 2010

 mris_inflate ../surf/lh.smoothwm ../surf/lh.inflated 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 ../surf/lh.inflated 


#-----------------------------------------
#@# Curvature Stats lh Thu Jul 15 16:58:44 EDT 2010

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/lh.curv.stats -F smoothwm FSL_MNI152_FS lh curv sulc 

#--------------------------------------------
#@# Sphere lh Thu Jul 15 16:58:52 EDT 2010

 mris_sphere -seed 1234 ../surf/lh.inflated ../surf/lh.sphere 

#--------------------------------------------
#@# Surf Reg lh Thu Jul 15 18:46:54 EDT 2010

 mris_register -curv ../surf/lh.sphere /usr/local/freesurfer/stable4//average/lh.average.curvature.filled.buckner40.tif ../surf/lh.sphere.reg 

#--------------------------------------------
#@# Jacobian white lh Thu Jul 15 21:14:03 EDT 2010

 mris_jacobian ../surf/lh.white ../surf/lh.sphere.reg ../surf/lh.jacobian_white 

#--------------------------------------------
#@# AvgCurv lh Thu Jul 15 21:14:07 EDT 2010

 mrisp_paint -a 5 /usr/local/freesurfer/stable4//average/lh.average.curvature.filled.buckner40.tif#6 ../surf/lh.sphere.reg ../surf/lh.avg_curv 

#-----------------------------------------
#@# Cortical Parc lh Thu Jul 15 21:14:10 EDT 2010

 mris_ca_label -aseg ../mri/aseg.mgz -seed 1234 FSL_MNI152_FS lh ../surf/lh.sphere.reg /usr/local/freesurfer/stable4//average/lh.curvature.buckner40.filled.desikan_killiany.2009-03-04.gcs ../label/lh.aparc.annot 

#-----------------------------------------
#@# Parcellation Stats lh Thu Jul 15 21:14:36 EDT 2010

 mris_anatomical_stats -mgz -f ../stats/lh.aparc.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab FSL_MNI152_FS lh 

#-----------------------------------------
#@# Cortical Parc 2 lh Thu Jul 15 21:14:45 EDT 2010

 mris_ca_label -aseg ../mri/aseg.mgz -seed 1234 FSL_MNI152_FS lh ../surf/lh.sphere.reg /usr/local/freesurfer/stable4//average/lh.destrieux.simple.2009-07-29.gcs ../label/lh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 lh Thu Jul 15 21:15:19 EDT 2010

 mris_anatomical_stats -mgz -f ../stats/lh.aparc.a2009s.stats -b -a ../label/lh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab FSL_MNI152_FS lh 

#--------------------------------------------
#@# Tessellate rh Thu Jul 15 21:15:30 EDT 2010

 mri_pretess ../mri/filled.mgz 127 ../mri/norm.mgz ../mri/filled-pretess127.mgz 


 mri_tessellate ../mri/filled-pretess127.mgz 127 ../surf/rh.orig.nofix 


 rm -f ../mri/filled-pretess127.mgz 


 mris_extract_main_component ../surf/rh.orig.nofix ../surf/rh.orig.nofix 

#--------------------------------------------
#@# Smooth1 rh Thu Jul 15 21:15:42 EDT 2010

 mris_smooth -nw -seed 1234 ../surf/rh.orig.nofix ../surf/rh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 rh Thu Jul 15 21:15:48 EDT 2010

 mris_inflate -no-save-sulc ../surf/rh.smoothwm.nofix ../surf/rh.inflated.nofix 

#--------------------------------------------
#@# QSphere rh Thu Jul 15 21:16:58 EDT 2010

 mris_sphere -q -seed 1234 ../surf/rh.inflated.nofix ../surf/rh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology rh Thu Jul 15 21:28:20 EDT 2010

 cp ../surf/rh.orig.nofix ../surf/rh.orig 


 cp ../surf/rh.inflated.nofix ../surf/rh.inflated 


 mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 FSL_MNI152_FS rh 


 mris_euler_number ../surf/rh.orig 


 mris_remove_intersection ../surf/rh.orig ../surf/rh.orig 


 rm ../surf/rh.inflated 

#--------------------------------------------
#@# Make Final Surf rh Thu Jul 15 22:00:00 EDT 2010

 mris_make_surfaces -noaparc -mgz -T1 brain.finalsurfs FSL_MNI152_FS rh 

#--------------------------------------------
#@# Surf Volume rh Thu Jul 15 23:03:37 EDT 2010

 mris_calc -o rh.area.mid rh.area add rh.area.pial 


 mris_calc -o rh.area.mid rh.area.mid div 2 


 mris_calc -o rh.volume rh.area.mid mul rh.thickness 

#--------------------------------------------
#@# Smooth2 rh Thu Jul 15 23:03:38 EDT 2010

 mris_smooth -n 3 -nw -seed 1234 ../surf/rh.white ../surf/rh.smoothwm 

#--------------------------------------------
#@# Inflation2 rh Thu Jul 15 23:03:44 EDT 2010

 mris_inflate ../surf/rh.smoothwm ../surf/rh.inflated 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 ../surf/rh.inflated 


#-----------------------------------------
#@# Curvature Stats rh Thu Jul 15 23:06:34 EDT 2010

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/rh.curv.stats -F smoothwm FSL_MNI152_FS rh curv sulc 

#--------------------------------------------
#@# Sphere rh Thu Jul 15 23:06:41 EDT 2010

 mris_sphere -seed 1234 ../surf/rh.inflated ../surf/rh.sphere 

#--------------------------------------------
#@# Surf Reg rh Fri Jul 16 01:03:51 EDT 2010

 mris_register -curv ../surf/rh.sphere /usr/local/freesurfer/stable4//average/rh.average.curvature.filled.buckner40.tif ../surf/rh.sphere.reg 

#--------------------------------------------
#@# Jacobian white rh Fri Jul 16 03:44:54 EDT 2010

 mris_jacobian ../surf/rh.white ../surf/rh.sphere.reg ../surf/rh.jacobian_white 

#--------------------------------------------
#@# AvgCurv rh Fri Jul 16 03:44:57 EDT 2010

 mrisp_paint -a 5 /usr/local/freesurfer/stable4//average/rh.average.curvature.filled.buckner40.tif#6 ../surf/rh.sphere.reg ../surf/rh.avg_curv 

#-----------------------------------------
#@# Cortical Parc rh Fri Jul 16 03:44:59 EDT 2010

 mris_ca_label -aseg ../mri/aseg.mgz -seed 1234 FSL_MNI152_FS rh ../surf/rh.sphere.reg /usr/local/freesurfer/stable4//average/rh.curvature.buckner40.filled.desikan_killiany.2009-03-04.gcs ../label/rh.aparc.annot 

#-----------------------------------------
#@# Parcellation Stats rh Fri Jul 16 03:45:27 EDT 2010

 mris_anatomical_stats -mgz -f ../stats/rh.aparc.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab FSL_MNI152_FS rh 

#-----------------------------------------
#@# Cortical Parc 2 rh Fri Jul 16 03:45:35 EDT 2010

 mris_ca_label -aseg ../mri/aseg.mgz -seed 1234 FSL_MNI152_FS rh ../surf/rh.sphere.reg /usr/local/freesurfer/stable4//average/rh.destrieux.simple.2009-07-29.gcs ../label/rh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 rh Fri Jul 16 03:46:11 EDT 2010

 mris_anatomical_stats -mgz -f ../stats/rh.aparc.a2009s.stats -b -a ../label/rh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab FSL_MNI152_FS rh 

#--------------------------------------------
#@# ASeg Stats Fri Jul 16 03:46:22 EDT 2010

 mri_segstats --seg mri/aseg.mgz --sum stats/aseg.stats --pv mri/norm.mgz --excludeid 0 --brain-vol-from-seg --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --surf-wm-vol --surf-ctx-vol --ctab /usr/local/freesurfer/stable4//ASegStatsLUT.txt --subject FSL_MNI152_FS 

#--------------------------------------------
#@# Cortical ribbon mask Fri Jul 16 03:57:13 EDT 2010

 mris_volmask --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon --save_distance FSL_MNI152_FS 

#-----------------------------------------
#@# AParc-to-ASeg Fri Jul 16 04:06:22 EDT 2010

 mri_aparc2aseg --s FSL_MNI152_FS --volmask 


 mri_aparc2aseg --s FSL_MNI152_FS --volmask --a2009s 

#-----------------------------------------
#@# WMParc Fri Jul 16 04:09:11 EDT 2010

 mri_aparc2aseg --s FSL_MNI152_FS --labelwm --hypo-as-wm --rip-unknown --volmask --o mri/wmparc.mgz --ctxseg aparc+aseg.mgz 


 mri_segstats --seg mri/wmparc.mgz --sum stats/wmparc.stats --pv mri/norm.mgz --excludeid 0 --brain-vol-from-seg --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --subject FSL_MNI152_FS --surf-wm-vol --ctab /usr/local/freesurfer/stable4//FreeSurferColorLUT.txt 

