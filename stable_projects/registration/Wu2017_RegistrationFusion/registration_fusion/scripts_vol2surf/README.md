## Registration Fusion (RF) Volume-to-Fsaverage Procedure

The RF approach (using GSP subjects) for volume-to-fsaverage mapping is implemented with the following steps.

Step 1: create index files of the volumetric atlas space

Step 2: project the index files to fsaverage (with RF-M3Z or RF-ANTs)

Step 3: compute average mapping across projected index files

----

## RF-ANTs Commands

Use the following commands in sequence to implement RF-ANTs approach for volume-to-surface mapping. For PBS scheduler users, you would need to wait for all the jobs submitted to be completed before calling the next command.

- for MNI152-to-fsaverage mapping:

    1) `../CBIG_RF_ANTsReg_prep.sh -p 'MNI152_orig'`

       or for PBS scheduler users (for example using circ-spool queue) `../CBIG_RF_ANTsReg_prep.sh -p 'MNI152_orig' -q circ-spool`

       This step can be skipped if you have previously run ANTs registration between the volumetric atlas and the individual subjects.
       
    2) `./CBIG_RF_step1_make_xyzIndex_volTemplate.sh -p 'MNI152_orig'`
    
    3) `./CBIG_RF_step2B_RF_ANTs_vol2fsaverage_proj.sh -p 'MNI152_orig'`
    
       or for PBS scheduler users (for example using circ-spool queue) `./CBIG_RF_step2B_RF_ANTs_vol2fsaverage_proj.sh -p 'MNI152_orig' -q circ-spool`

       if you want to use your own warps (meaning you have skipped step 1), you can pass in the directory containing your own warps with the `-w` option.
    
    4) `./CBIG_RF_step3_compute_vol2fsaverage_avgMapping.sh -p 'MNI152_orig' -c 0`

- for Colin27-to-fsaverage mapping:

    1) `../CBIG_RF_ANTsReg_prep.sh -p 'Colin27_orig'`
    
       or for PBS scheduler users (for example using circ-spool queue) `../CBIG_RF_ANTsReg_prep.sh -p 'Colin27_orig' -q circ-spool`

       This step can be skipped if you have previously run ANTs registration between the volumetric atlas and
the individual subjects.
       
    2) `./CBIG_RF_step1_make_xyzIndex_volTemplate.sh -p 'Colin27_orig' -o $(pwd)/results/index_Colin27`
    
    3) `./CBIG_RF_step2B_RF_ANTs_vol2fsaverage_proj.sh -p 'Colin27_orig' -i $(pwd)/results/index_Colin27`
    
       or for PBS scheduler users (for example using circ-spool queue) `./CBIG_RF_step2B_RF_ANTs_vol2fsaverage_proj.sh -p 'Colin27_orig' -i $(pwd)/results/index_Colin27 -q circ-spool`
       if you want to use your own warps (meaning you have skipped step 1), you can pass in the directory containing your own warps with the `-w` option.
    
    4) `./CBIG_RF_step3_compute_vol2fsaverage_avgMapping.sh -p 'Colin27_orig' -i $(pwd)/results/index_fsaverage -c 0`

----

## RF-M3Z Commands

Use the following commands in sequence to implement RF-M3Z approach for volume-to-surface mapping. For PBS scheduler users, you would need to wait for all the jobs submitted to be completed before calling the next command.

- for MNI152-to-fsaverage mapping:

    1) `./CBIG_RF_step1_make_xyzIndex_volTemplate.sh -p 'MNI152_norm'`
    
    2) `./CBIG_RF_step2A_RF_M3Z_vol2fsaverage_proj.sh -p 'MNI152_norm'` 
    
       or for PBS scheduler users (for example using circ-spool queue) `./CBIG_RF_step2A_RF_M3Z_vol2fsaverage_proj.sh -p 'MNI152_norm' -q circ-spool`
       
    3) `./CBIG_RF_step3_compute_vol2fsaverage_avgMapping.sh -p 'MNI152_norm' -r 'RF_M3Z' -c 0`
    
- for Colin27-to-fsaverage mapping:

    1) `./CBIG_RF_step1_make_xyzIndex_volTemplate.sh -p 'Colin27_norm' -o $(pwd)/results/index_Colin27`
    
    2) `./CBIG_RF_step2A_RF_M3Z_vol2fsaverage_proj.sh -p 'Colin27_norm' -i $(pwd)/results/index_Colin27 -s 'SPM_Colin27_FS4.5.0'` 
    
       or for PBS scheduler users (for example using circ-spool queue) `./CBIG_RF_step2A_RF_M3Z_vol2fsaverage_proj.sh -p 'Colin27_norm' -i $(pwd)/results/index_Colin27 -s 'SPM_Colin27_FS4.5.0' -q circ-spool`
       
    3) `./CBIG_RF_step3_compute_vol2fsaverage_avgMapping.sh -p 'Colin27_norm' -i $(pwd)/results/index_fsaverage -r 'RF_M3Z' -c 0`



    
    
