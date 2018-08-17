This README includes the instructions on how to run CBIG_gwMRF unit test. Notice that all filenames and directories in this unit test **work for CBIG lab only**.

----

References
==========
+ Schaefer A, Kong R, Gordon EM, Laumann TO, Zuo XN, Holmes AJ, Eickhoff SB, Yeo BTT. [**Local-Global parcellation of the human cerebral cortex from intrinsic functional connectivity MRI**](http://people.csail.mit.edu/ythomas/publications/2018LocalGlobal-CerebCor.pdf), *Cerebral Cortex*, 29:3095-3114, 2018

----

Data
====
We use 10 preprocessed subjects from the GSP dataset for this unit test. The data can be found under `/mnt/eql/yeo1/CBIG_private_unit_tests_data/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/data`.

----

Run
===
In terminal, call the unit test script:
```
$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/unit_tests/scripts/CBIG_gwMRF_unit_test.sh <your_output_folder>
```

This will submit a job to `circ-spool` to create two parcellation randomizations with 50 left and 50 right hemisphere cluster in your output folder.
 
Then the script will automatically call function `CBIG_gwMRF_check_unit_test_result.m` to compare your results with the results in our `ref_dir` (/mnt/eql/yeo1/CBIG_private_unit_tests_data/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/results).

A log file will be created in `<your_output_folder>/logs`, containing detailed information about whether the unit test is `[PASSED]` or `[FAILED]`. 

----

Input
=====
You need to specify `<your_output_folder>` when you call `CBIG_gwMRF_unit_test.sh`. Unit test parcellation results and log files will be stored there.

----

Output
======
- Two concatenated time matrices in `<your_output_folder>/time_data`.
 
- Two premultiplied product matrices in `<your_output_folder>/mult_mat`.

- Two sets of random parcellation results in `<your_output_folder>/clustering`. Each randomization contains two mat files.
  For example for seed 1:
  1) *_seed_1_Energy.mat
  2) *_seed_1.mat

- A log file in `<your_output_folder>/logs`. This log file contains the comparison results between your above mentioned output (time matrices, product matrices, *_seed_*.mats) and the output in our `ref_dir`.
  
  *[IMPORTANT]*

  All you need to do is to check this log file after the unit test is ended. A `[FAILED]` message indicates something wrong with the codes. The unit test is successful only if all messages are `[PASSED]`.

----

Bugs and Questions
==================
Please contact Alexander Schaefer at alexschaefer83@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
