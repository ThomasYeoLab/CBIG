# Brain Extraction with FSL BET -- General

**Note** This is the general procedure, which is easier to be applied to a new dataset. To replicate our PNAS results or view a simplified version for concept/procedure understanding, please refer to folder `../replicatePNAS/` instead.

---
## Master Script `CBIG_runBET.sh`

This script runs all BET jobs serially or submits them to a job queue if you have a cluster.

Before running this script, you should

* Prepare a text file whose lines are paths to NIfTI volumes that you want to perform BET on
* Prepare another text file specifying output names for *all* volumes, if you want to use names different from the original filenames. This is useful when you want your BET outputs to have customized filenames (e.g., `0002_sc`) rather than build on top of the original filenames (e.g., `ADNI_011_S_0002_MR_MPR__GradWarp__B1_Correction__N3__Scaled_Br_20070108225928642_S9107_I35475`).
* Customize BET parameters by editing the job scripts (`BETParam*.sh`)

`../replicatePNAS/CBIG_runBET_wrapper.sh` demonstrates how to call `CBIG_runBET.sh`. You may also do `./CBIG_runBET.sh` to see detailed usage.

----
## Job Scripts `CBIG_BETParam*.sh`

These are job scripts (to be run or submitted by `CBIG_runBET.sh`), where you can adapt the BET procedure or just parameters based on deficiencies or errors from the previous round of BET.

More specifically, `CBIG_BETParam1.sh` specifies the procedure or parameters for the initial round of BET, which is performed on all images. `CBIG_BETParam2.sh` specifies a different procedure or different parameters for the second round of BET, which is performed only on the images that the first round has failed on. You may need to do several rounds of BET until the BET results of all images are satisfactory.

----
## Quality Inspection Script `CBIG_QC_slices.sh`

This script overlays each extracted brain on top of the corresponding original volume so that you can inspect the BET results and decide which results are unsatisfactory and should go for another round of BET. Open `index.html` in a web browser to view the results.

To see detailed usage, do `./CBIG_QC_slices.sh`.

----
## Helper Script for Quality Control `CBIG_QC_getOrigPaths.sh`

In the quality inspection, you will note down the images that need another round of BET with different parameters. During this procedure, it is much more convinient to note down the brain filenames (e.g., `0002_sc_brain`) from `index.html` than to look up the original filenames (e.g., `ADNI_011_S_0002_MR_MPR__GradWarp__B1_Correction__N3__Scaled_Br_20070108225928642_S9107_I35475.nii`), although we do need the original paths for the next round of BET. 

This helper script takes in a text file whose lines are brain filenames and outputs a new text file whose lines are orginal paths. This new text file can then be fed into `CBIG_runBET` for the next round of BET.

To see detailed usage, do `./CBIG_QC_getOrigPaths.sh`. 

----
## Summarization Function `CBIG_selectBETResults.m`

When you are satisfied with all BET results, call this MATLAB function to generate a list of final BET results. For each image, the function selects the best BET result, which is essentially the result in the latest round of BET. For instance, if an image is BET'ed in rounds 1, 2 and 3 (implying you sent it to Round 3 because you were unsatisfied with the Round 1 and 2 results!), then the function takes the round 3 result as the final BET result for this image.

The output file (`brainList.txt`) is in the same folder as the BET output folders and should be used for the subsequent steps (VBM).

To see detailed usage, do `help CBIG_selectBETResults` in MATLAB.

----
## Overall Procedure

Here is an example workflow:

1. First round of BET on all images
	* Call `CBIG_runBET.sh` with `CBIG_BETParam1.sh` as one of the arguments
2. Quality check of the first round
	1. Run `CBIG_QC_slices.sh` to screenshot brains overlaid on top of the original images
	2. Copy and paste filenames of unsatisfactory brains (e.g., 0002_sc_brain) into a text file
	3. Run `CBIG_QC_getOrigPaths.sh` over the text file to translate brain filenames into paths to the original images
3. Second round of BET on subset from the previous round
	1. Edit `CBIG_BETParam2.sh` based on the manifested deficiency or errors of `CBIG_BETParam1.sh`
	2. Call `CBIG_runBET.sh` with `CBIG_BETParam2.sh` as one of the arguments
4. Quality check of the second round
5. Repeat Steps 3 and 4 until all BET results are satisfactory
6. Call `CBIG_selectBETResults.m` in MATLAB to get a list of final BET results for the subsequent VBM steps

