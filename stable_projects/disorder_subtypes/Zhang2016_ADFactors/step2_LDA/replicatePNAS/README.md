# LDA -- Specific to Our PNAS Paper

**NOTE** If you encounter `Segmentation fault` at the LDA line, it is often the case that your argument strings are too long. Please shorten them by, for example, placing the input files together with the script (so that you can reference them with a very short string, i.e., `./`).

This folder allows you to replicate<sup>1</sup> our LDA results reported in our PNAS paper.

The wrapper first converts VBM results from the previous step into "brain documents" that can be analyzed by LDA, next estimates the LDA model parameters for K = 2, 3 and 4 with the 188 ADNI-1 AD patients, then visualizes the extracted atrophy factors as probabilistic atrophy maps, and finally infers factor compositions of participants other than the 188 used to estimate the model parameters.

Folder `inputs_brain2doc/` provides the input files that we used during brain-to-document conversion. They can be used as format references, which you should follow while generating your own input files. 

<sup>1</sup>Not exact replication, because results vary across different random seeds. As such, we have provided our LDA outputs (used in the PNAS paper) in `../../ADNIDataRelease/inferNew/files_LDA/`. However, in our experience, the results are highly similar despite different random seeds. 

----
## How do lines in the factor composition output file correspond to the images?

The principle here is that lines in the factor composition output file always match with order of the input documents (images). Hence, if you input a filename list to `../lib/CBIG_brain2doc.m`, lines in the factor composition output file will correspond to that list. That is, the first line of factor compositions will correspond to the first image listed in the filename list. Otherwise, they will match the order of the `concatOrder` file.
