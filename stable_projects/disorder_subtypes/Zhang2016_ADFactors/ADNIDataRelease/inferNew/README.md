## Looking for the probabilistic atrophy maps?

Please go to folders `files_LDA/model_K2`, `files_LDA/model_K3` and `files_LDA/model_K4` for the two-, three- and four-factor models, respectively. You can also view or download all of these atrophy maps online at [http://neurovault.org/collections/1917/](http://neurovault.org/collections/1917/).

----

# Inferring Factor Compositions of Unseen Subjects

**NOTE** If you encounter `Segmentation fault` at the LDA line, it is often the case that your argument strings are too long. Please shorten them by, for example, placing the input files together with the script (so that you can reference them with a very short string, i.e., `./`).

---
## Wrapper `CBIG_inferNew_wrapper.sh`

This is a wrapper for `CBIG_inferNew.sh` aiming to demonstrate how to call `CBIG_inferNew.sh`.

----
## Main Script `CBIG_inferNew.sh`

It essentially calls the functions in `../../step2_LDA/lib/` to perform VBM on the unseen (new) subjects, then convert the VBM results into "brain documents", and finally infer factor compositions of these unseen subjects.

`CBIG_reorderNuisance.m` is simply a helper function that reorders lines in the nuisance variable file. It will be automatically called in `CBIG_inferNew.sh`, so you do not need to worry about it. The following explanation is only for those interested in knowing what it does exactly and why it is necessary. Lines in the nuisane variable file match with the list of input scans, but not necessarily with the concatenation order of the VBM results. Hence, this helper function reorders the nuisance variable lines so that they match with order of the VBM results, which, of course, is a must for correct regression.

----
## Necessary Files

Besides functions and scripts, some necessary VBM and LDA files are also released in folders `files_VBM` and `files_LDA`, respectively.

----
## How do lines in the factor composition output file correspond to the images?

The principle here is that lines in the factor composition output file always match with order of the input documents (images). Hence, if you input a filename list to `../lib/CBIG_brain2doc.m`, lines in the factor composition output file will correspond to that list. That is, the first line of factor compositions will correspond to the first image listed in the filename list. Otherwise, they will match the order of the `concatOrder` file.
