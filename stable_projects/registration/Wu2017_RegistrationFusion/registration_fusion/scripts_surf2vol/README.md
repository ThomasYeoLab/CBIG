## Registration Fusion (RF) Fsaverage-to-Volume Procedure

The RF approach (using GSP subjects) for fsaverage-to-volume mapping is implemented with the following steps.

Step 1: create index files of fsaverage surface (which are projected to each GSP subject's surface space)

Step 2: project the index files to the volumetric atlas (with RF-M3Z or RF-ANTs)

Step 3: compute average mapping across projected index files & generate corresponding liberal cortex mask

----

## RF-ANTs Commands

Use the following commands in sequence to implement RF-ANTs approach for surface-to-volume mapping. For PBS scheduler users, you would need to wait for all the jobs submitted to be completed before calling the next command.

- for fsaverage-to-MNI152 mapping:

    1) `../CBIG_RF_ANTsReg_prep.sh -p 'MNI152_orig'`

       or for PBS scheduler users (for example using circ-spool queue) `../CBIG_RF_ANTsReg_prep.sh -p 'MNI152_orig' -q circ-spool`

       This step can be skipped if you have previously run ANTs registration between the volumetric atlas and
the individual subjects.
       
    2) `./CBIG_RF_step1_make_xyzIndex_fsaverage.sh -n 0`
    
    3) `./CBIG_RF_step2B_RF_ANTs_fsaverage2vol_proj.sh -p 'MNI152_orig' -s 'FSL_MNI152_FS4.5.0'`
    
       or for PBS scheduler users (for example using circ-spool queue) `./CBIG_RF_step2B_RF_ANTs_fsaverage2vol_proj.sh -p 'MNI152_orig' -s 'FSL_MNI152_FS4.5.0'`
       if you want to use your own warps (meaning you have skipped step 1), you can pass in the directory containing your own warps with the `-w` option.
    
    4) `./CBIG_RF_step3_compute_fsaverage2vol_avgMapping.sh -s 'FSL_MNI152_FS4.5.0' -c 0`

- for fsaverage-to-Colin27 mapping:

    1) `../CBIG_RF_ANTsReg_prep.sh -p 'Colin27_orig'`
    
       or for PBS scheduler users (for example using circ-spool queue) `../CBIG_RF_ANTsReg_prep.sh -p 'Colin27_orig' -q circ-spool`

       This step can be skipped if you have previously run ANTs registration between the volumetric atlas and
the individual subjects.
       
    2) `./CBIG_RF_step1_make_xyzIndex_fsaverage.sh -n 0`
    
    3) `./CBIG_RF_step2B_RF_ANTs_fsaverage2vol_proj.sh -p 'Colin27_orig' -s 'SPM_Colin27_FS4.5.0'`
    
       or for PBS scheduler users (for example using circ-spool queue) `./CBIG_RF_step2B_RF_ANTs_fsaverage2vol_proj.sh -p 'Colin27_orig' -s 'SPM_Colin27_FS4.5.0' -q circ-spool`
       if you want to use your own warps (meaning you have skipped step 1), you can pass in the directory containing your own warps with the `-w` option.
    
    4) `./CBIG_RF_step3_compute_fsaverage2vol_avgMapping.sh -s 'SPM_Colin27_FS4.5.0' -c 0`

----

## RF-M3Z Commands

Use the following commands in sequence to implement RF-M3Z approach for surface-to-volume mapping. For PBS scheduler users, you would need to wait for all the jobs submitted to be completed before calling the next command.

- for fsaverage-to-MNI152 mapping:

    1) `./CBIG_RF_step1_make_xyzIndex_fsaverage.sh -n 0`
    
    2) `./CBIG_RF_step2A_RF_M3Z_fsaverage2vol_proj.sh -s 'FSL_MNI152_FS4.5.0'` 
    
       or for PBS scheduler users (for example using circ-spool queue) `./CBIG_RF_step2A_RF_M3Z_fsaverage2vol_proj.sh -s 'FSL_MNI152_FS4.5.0' -q circ-spool`
       
    3) `./CBIG_RF_step3_compute_fsaverage2vol_avgMapping.sh -s 'FSL_MNI152_FS4.5.0' -r 'RF_M3Z' -c 0`
    
- for fsaverage-to-Colin27 mapping:

    1) `./CBIG_RF_step1_make_xyzIndex_fsaverage.sh -n 0`
    
    2) `./CBIG_RF_step2A_RF_M3Z_fsaverage2vol_proj.sh -s 'SPM_Colin27_FS4.5.0'` 
    
       or for PBS scheduler users (for example using circ-spool queue) `./CBIG_RF_step2A_RF_M3Z_fsaverage2vol_proj.sh -s 'SPM_Colin27_FS4.5.0' -q circ-spool`
       
    3) `./CBIG_RF_step3_compute_fsaverage2vol_avgMapping.sh -s 'SPM_Colin27_FS4.5.0' -r 'RF_M3Z' -c 0`
    

    
    
