
 mri_convert /mnt/eql/p1/users/external/jianxiaow/projects/Git_codes/data/templates/colin27T1_seg.nii /mnt/eql/p1/users/external/jianxiaow/projects/Git_codes/data/Colin27_45/mri/orig/001.mgz 

#--------------------------------------------
#@# MotionCor Tue Jun 21 18:20:41 SGT 2016

 cp /mnt/eql/p1/users/external/jianxiaow/projects/Git_codes/data/Colin27_45/mri/orig/001.mgz /mnt/eql/p1/users/external/jianxiaow/projects/Git_codes/data/Colin27_45/mri/rawavg.mgz 


 mri_convert /mnt/eql/p1/users/external/jianxiaow/projects/Git_codes/data/Colin27_45/mri/rawavg.mgz /mnt/eql/p1/users/external/jianxiaow/projects/Git_codes/data/Colin27_45/mri/orig.mgz --conform 


 mri_add_xform_to_header -c /mnt/eql/p1/users/external/jianxiaow/projects/Git_codes/data/Colin27_45/mri/transforms/talairach.xfm /mnt/eql/p1/users/external/jianxiaow/projects/Git_codes/data/Colin27_45/mri/orig.mgz /mnt/eql/p1/users/external/jianxiaow/projects/Git_codes/data/Colin27_45/mri/orig.mgz 

#--------------------------------------------
#@# Nu Intensity Correction Tue Jun 21 18:20:48 SGT 2016

 mri_nu_correct.mni --i orig.mgz --o nu.mgz --n 2 

#--------------------------------------------
#@# Talairach Tue Jun 21 18:21:11 SGT 2016

 talairach_avi --i nu.mgz --xfm transforms/talairach.auto.xfm 

#--------------------------------------------
#@# Talairach Failure Detection Tue Jun 21 18:21:46 SGT 2016

 talairach_afd -T 0.005 -xfm transforms/talairach.xfm 


 awk -f /apps/arch/Linux_x86_64/freesurfer/4.5.0/bin/extract_talairach_avi_QA.awk /mnt/eql/p1/users/external/jianxiaow/projects/Git_codes/data/Colin27_45/mri/transforms/talairach_avi.log 

#--------------------------------------------
#@# Intensity Normalization Tue Jun 21 18:21:46 SGT 2016

 mri_normalize -g 1 nu.mgz T1.mgz 

#--------------------------------------------
#@# Skull Stripping Tue Jun 21 18:23:52 SGT 2016

 mri_watershed -keep brainmask.auto.mgz brainmask.mgz brainmask.mgz -T1 -brain_atlas /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta T1.mgz brainmask.auto.mgz 


INFO: brainmask.mgz already exists!
The new brainmask.auto.mgz will not be copied to brainmask.mgz.
This is done to retain any edits made to brainmask.mgz.
Add the -clean-bm flag to recon-all to overwrite brainmask.mgz.

#-------------------------------------
#@# EM Registration Tue Jun 21 18:24:48 SGT 2016

 mri_em_register -mask brainmask.mgz nu.mgz /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/RB_all_2008-03-26.gca transforms/talairach.lta 

#--------------------------------------
#@# CA Normalize Tue Jun 21 18:41:24 SGT 2016

 mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/RB_all_2008-03-26.gca transforms/talairach.lta norm.mgz 

#--------------------------------------
#@# CA Reg Tue Jun 21 18:43:00 SGT 2016

 mri_ca_register -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/RB_all_2008-03-26.gca transforms/talairach.m3z 

#--------------------------------------
#@# CA Reg Inv Tue Jun 21 21:27:32 SGT 2016

 mri_ca_register -invert-and-save transforms/talairach.m3z 

#--------------------------------------
#@# Remove Neck Tue Jun 21 21:28:14 SGT 2016

 mri_remove_neck -radius 25 nu.mgz transforms/talairach.m3z /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/RB_all_2008-03-26.gca nu_noneck.mgz 

#--------------------------------------
#@# SkullLTA Tue Jun 21 21:29:06 SGT 2016

 mri_em_register -skull -t transforms/talairach.lta nu_noneck.mgz /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta 

#--------------------------------------
#@# SubCort Seg Tue Jun 21 22:02:02 SGT 2016

 mri_seg_diff --seg1 aseg.auto.mgz --seg2 aseg.mgz --diff aseg.manedit.mgz 


 mri_ca_label -align -nobigventricles norm.mgz transforms/talairach.m3z /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/RB_all_2008-03-26.gca aseg.auto_noCCseg.mgz 


 mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz Colin27_45 

#--------------------------------------
#@# Merge ASeg Tue Jun 21 22:20:33 SGT 2016

 cp aseg.auto.mgz aseg.mgz 

#--------------------------------------------
#@# Intensity Normalization2 Tue Jun 21 22:20:33 SGT 2016

 mri_normalize -aseg aseg.mgz -mask brainmask.mgz norm.mgz brain.mgz 

#--------------------------------------------
#@# Mask BFS Tue Jun 21 22:23:26 SGT 2016

 mri_mask -T 5 brain.mgz brainmask.mgz brain.finalsurfs.mgz 

#--------------------------------------------
#@# WM Segmentation Tue Jun 21 22:23:28 SGT 2016

 cp wm.mgz wm.seg.mgz 


 mri_segment -keep brain.mgz wm.seg.mgz 


 mri_edit_wm_with_aseg -keep-in wm.seg.mgz brain.mgz aseg.mgz wm.asegedit.mgz 


 mri_pretess -keep wm.asegedit.mgz wm norm.mgz wm.mgz 

#--------------------------------------------
#@# Fill Tue Jun 21 22:25:15 SGT 2016

 mri_fill -a ../scripts/ponscc.cut.log -xform transforms/talairach.lta -segmentation aseg.auto_noCCseg.mgz wm.mgz filled.mgz 

#--------------------------------------------
#@# Tessellate lh Tue Jun 21 22:26:03 SGT 2016

 mri_pretess ../mri/filled.mgz 255 ../mri/norm.mgz ../mri/filled-pretess255.mgz 


 mri_tessellate ../mri/filled-pretess255.mgz 255 ../surf/lh.orig.nofix 


 rm -f ../mri/filled-pretess255.mgz 


 mris_extract_main_component ../surf/lh.orig.nofix ../surf/lh.orig.nofix 

#--------------------------------------------
#@# Smooth1 lh Tue Jun 21 22:26:10 SGT 2016

 mris_smooth -nw -seed 1234 ../surf/lh.orig.nofix ../surf/lh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 lh Tue Jun 21 22:26:14 SGT 2016

 mris_inflate -no-save-sulc ../surf/lh.smoothwm.nofix ../surf/lh.inflated.nofix 

#--------------------------------------------
#@# QSphere lh Tue Jun 21 22:27:04 SGT 2016

 mris_sphere -q -seed 1234 ../surf/lh.inflated.nofix ../surf/lh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology lh Tue Jun 21 22:34:52 SGT 2016

 cp ../surf/lh.orig.nofix ../surf/lh.orig 


 cp ../surf/lh.inflated.nofix ../surf/lh.inflated 


 mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 Colin27_45 lh 


 mris_euler_number ../surf/lh.orig 


 mris_remove_intersection ../surf/lh.orig ../surf/lh.orig 


 rm ../surf/lh.inflated 

#--------------------------------------------
#@# Make Final Surf lh Tue Jun 21 23:05:05 SGT 2016

 mris_make_surfaces -noaparc -mgz -T1 brain.finalsurfs Colin27_45 lh 

#--------------------------------------------
#@# Surf Volume lh Wed Jun 22 00:19:01 SGT 2016

 mris_calc -o lh.area.mid lh.area add lh.area.pial 


 mris_calc -o lh.area.mid lh.area.mid div 2 


 mris_calc -o lh.volume lh.area.mid mul lh.thickness 

#--------------------------------------------
#@# Smooth2 lh Wed Jun 22 00:19:01 SGT 2016

 mris_smooth -n 3 -nw -seed 1234 ../surf/lh.white ../surf/lh.smoothwm 

#--------------------------------------------
#@# Inflation2 lh Wed Jun 22 00:19:05 SGT 2016

 mris_inflate ../surf/lh.smoothwm ../surf/lh.inflated 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 ../surf/lh.inflated 


#-----------------------------------------
#@# Curvature Stats lh Wed Jun 22 00:21:28 SGT 2016

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/lh.curv.stats -F smoothwm Colin27_45 lh curv sulc 

#--------------------------------------------
#@# Sphere lh Wed Jun 22 00:21:33 SGT 2016

 mris_sphere -seed 1234 ../surf/lh.inflated ../surf/lh.sphere 

#--------------------------------------------
#@# Surf Reg lh Wed Jun 22 01:38:20 SGT 2016

 mris_register -curv ../surf/lh.sphere /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/lh.average.curvature.filled.buckner40.tif ../surf/lh.sphere.reg 

#--------------------------------------------
#@# Jacobian white lh Wed Jun 22 03:47:46 SGT 2016

 mris_jacobian ../surf/lh.white ../surf/lh.sphere.reg ../surf/lh.jacobian_white 

#--------------------------------------------
#@# AvgCurv lh Wed Jun 22 03:47:49 SGT 2016

 mrisp_paint -a 5 /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/lh.average.curvature.filled.buckner40.tif#6 ../surf/lh.sphere.reg ../surf/lh.avg_curv 

#-----------------------------------------
#@# Cortical Parc lh Wed Jun 22 03:47:51 SGT 2016

 mris_ca_label -aseg ../mri/aseg.mgz -seed 1234 Colin27_45 lh ../surf/lh.sphere.reg /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/lh.curvature.buckner40.filled.desikan_killiany.2009-03-04.gcs ../label/lh.aparc.annot 

#-----------------------------------------
#@# Parcellation Stats lh Wed Jun 22 03:48:16 SGT 2016

 mris_anatomical_stats -mgz -f ../stats/lh.aparc.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab Colin27_45 lh 

#-----------------------------------------
#@# Cortical Parc 2 lh Wed Jun 22 03:48:22 SGT 2016

 mris_ca_label -aseg ../mri/aseg.mgz -seed 1234 Colin27_45 lh ../surf/lh.sphere.reg /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/lh.destrieux.simple.2009-07-29.gcs ../label/lh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 lh Wed Jun 22 03:48:58 SGT 2016

 mris_anatomical_stats -mgz -f ../stats/lh.aparc.a2009s.stats -b -a ../label/lh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab Colin27_45 lh 

#--------------------------------------------
#@# Tessellate rh Wed Jun 22 03:49:04 SGT 2016

 mri_pretess ../mri/filled.mgz 127 ../mri/norm.mgz ../mri/filled-pretess127.mgz 


 mri_tessellate ../mri/filled-pretess127.mgz 127 ../surf/rh.orig.nofix 


 rm -f ../mri/filled-pretess127.mgz 


 mris_extract_main_component ../surf/rh.orig.nofix ../surf/rh.orig.nofix 

#--------------------------------------------
#@# Smooth1 rh Wed Jun 22 03:49:12 SGT 2016

 mris_smooth -nw -seed 1234 ../surf/rh.orig.nofix ../surf/rh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 rh Wed Jun 22 03:49:17 SGT 2016

 mris_inflate -no-save-sulc ../surf/rh.smoothwm.nofix ../surf/rh.inflated.nofix 

#--------------------------------------------
#@# QSphere rh Wed Jun 22 03:50:07 SGT 2016

 mris_sphere -q -seed 1234 ../surf/rh.inflated.nofix ../surf/rh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology rh Wed Jun 22 03:58:15 SGT 2016

 cp ../surf/rh.orig.nofix ../surf/rh.orig 


 cp ../surf/rh.inflated.nofix ../surf/rh.inflated 


 mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 Colin27_45 rh 


 mris_euler_number ../surf/rh.orig 


 mris_remove_intersection ../surf/rh.orig ../surf/rh.orig 


 rm ../surf/rh.inflated 

#--------------------------------------------
#@# Make Final Surf rh Wed Jun 22 04:10:27 SGT 2016

 mris_make_surfaces -noaparc -mgz -T1 brain.finalsurfs Colin27_45 rh 

#--------------------------------------------
#@# Surf Volume rh Wed Jun 22 05:22:24 SGT 2016

 mris_calc -o rh.area.mid rh.area add rh.area.pial 


 mris_calc -o rh.area.mid rh.area.mid div 2 


 mris_calc -o rh.volume rh.area.mid mul rh.thickness 

#--------------------------------------------
#@# Smooth2 rh Wed Jun 22 05:22:24 SGT 2016

 mris_smooth -n 3 -nw -seed 1234 ../surf/rh.white ../surf/rh.smoothwm 

#--------------------------------------------
#@# Inflation2 rh Wed Jun 22 05:22:29 SGT 2016

 mris_inflate ../surf/rh.smoothwm ../surf/rh.inflated 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 ../surf/rh.inflated 


#-----------------------------------------
#@# Curvature Stats rh Wed Jun 22 05:24:51 SGT 2016

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/rh.curv.stats -F smoothwm Colin27_45 rh curv sulc 

#--------------------------------------------
#@# Sphere rh Wed Jun 22 05:24:56 SGT 2016

 mris_sphere -seed 1234 ../surf/rh.inflated ../surf/rh.sphere 

#--------------------------------------------
#@# Surf Reg rh Wed Jun 22 06:55:10 SGT 2016

 mris_register -curv ../surf/rh.sphere /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/rh.average.curvature.filled.buckner40.tif ../surf/rh.sphere.reg 

#--------------------------------------------
#@# Jacobian white rh Wed Jun 22 09:12:58 SGT 2016

 mris_jacobian ../surf/rh.white ../surf/rh.sphere.reg ../surf/rh.jacobian_white 

#--------------------------------------------
#@# AvgCurv rh Wed Jun 22 09:13:00 SGT 2016

 mrisp_paint -a 5 /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/rh.average.curvature.filled.buckner40.tif#6 ../surf/rh.sphere.reg ../surf/rh.avg_curv 

#-----------------------------------------
#@# Cortical Parc rh Wed Jun 22 09:13:02 SGT 2016

 mris_ca_label -aseg ../mri/aseg.mgz -seed 1234 Colin27_45 rh ../surf/rh.sphere.reg /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/rh.curvature.buckner40.filled.desikan_killiany.2009-03-04.gcs ../label/rh.aparc.annot 

#-----------------------------------------
#@# Parcellation Stats rh Wed Jun 22 09:13:27 SGT 2016

 mris_anatomical_stats -mgz -f ../stats/rh.aparc.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab Colin27_45 rh 

#-----------------------------------------
#@# Cortical Parc 2 rh Wed Jun 22 09:13:33 SGT 2016

 mris_ca_label -aseg ../mri/aseg.mgz -seed 1234 Colin27_45 rh ../surf/rh.sphere.reg /apps/arch/Linux_x86_64/freesurfer/4.5.0/average/rh.destrieux.simple.2009-07-29.gcs ../label/rh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 rh Wed Jun 22 09:14:07 SGT 2016

 mris_anatomical_stats -mgz -f ../stats/rh.aparc.a2009s.stats -b -a ../label/rh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab Colin27_45 rh 

#--------------------------------------------
#@# ASeg Stats Wed Jun 22 09:14:14 SGT 2016

 mri_segstats --seg mri/aseg.mgz --sum stats/aseg.stats --pv mri/norm.mgz --excludeid 0 --brain-vol-from-seg --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --subject Colin27_45 --surf-wm-vol --ctab /apps/arch/Linux_x86_64/freesurfer/4.5.0/ASegStatsLUT.txt 

#--------------------------------------------
#@# Cortical ribbon mask Wed Jun 22 09:22:01 SGT 2016

 mris_volmask --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon --save_distance Colin27_45 

#-----------------------------------------
#@# AParc-to-ASeg Wed Jun 22 09:40:40 SGT 2016

 mri_aparc2aseg --s Colin27_45 --volmask 


 mri_aparc2aseg --s Colin27_45 --volmask --a2009s 

#-----------------------------------------
#@# WMParc Wed Jun 22 09:43:54 SGT 2016

 mri_aparc2aseg --s Colin27_45 --labelwm --hypo-as-wm --rip-unknown --volmask --o mri/wmparc.mgz --ctxseg aparc+aseg.mgz 


 mri_segstats --seg mri/wmparc.mgz --sum stats/wmparc.stats --pv mri/norm.mgz --excludeid 0 --brain-vol-from-seg --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --subject Colin27_45 --surf-wm-vol --ctab /apps/arch/Linux_x86_64/freesurfer/4.5.0/FreeSurferColorLUT.txt 

