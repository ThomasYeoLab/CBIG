## Background
* Head motion in the fMRI scanner corrupts fMRI data
* To minimize the quality degradation due to head motion, motion censoring (removing high motion frames from fMRI data) is commonly used
* In motion censoring, we typically compute the framewise motion estimates and use the motion estimates to identify high motion frames
* However, respiration could inflate the motion estimates, thus resulting in unwanted censoring of good-quality frames
* Our processing pipeline allows users to remove the respiratory pseudomotion from motion estimates (respiratory pseudomotion filtering)

## What is respiratory pseudomotion (Power et.al. 2019)
* Chest movement during respiration generate small perturbations of main magnetic field(B0) of scanner
* This perturbation of B0 causes shift in reconstructed image in the phase encoding direction
* This problem gets worse in multiband scanners for their higher sampling frequencies
* This reconstruction shift is different from true motion, as in it doesn't cause imaging quality issues like true motion
* Thus, we would like to filter out respiratory pseudomotion from motion estimates to save more data

## Filtering of respiratory pseudomotion (multi-band)
* For multi-band data, we just apply a bandstop filtering on the motion estimates (Fair et.al. 2020). For example, if we know the respiratory frequency lies between 0.31 to 0.43 Hz, we can set the stopband to be 0.31 to 0.43 Hz
* The stopband should vary for different datasets because different population could have different respiration frequency. For example, children have higher respiration rates than adults
* Ideally, the stopband should be calculated using respiration data. For example, Fair et.al. 2020 used respiration data of ABCD participants, plotted the histogram of the respiration frequency, and defined the stopband as between 2nd and 3rd quantiles of the respiration frequency histogram
* Alternatively, the stopband could be inferred from the age of participants given that respiration rates varies with age. Fleming et.al. 2011 listed the respiration rates for children of different age groups in Supplementary Table 4

## Filtering of respiratory pseudomotion (single-band)
* Single-band data usually have lower sampling frequency (i.e., higher TR). As a result, the Nyquist frequency could be lower than respiration frequency and the respiration frequency motion aliases into other frequency bands
* Gratton et.al. 2020 proposed to use a low-pass filter and remove motion frequency > 0.1 Hz (high frequency motion)

## Implementation in CBIG preprocessing pipeline
* In CBIG preprocessing pipeline, we have options to add respiratory pseudomotion filtering during the motion correction
* To perform respiratory pseudomotion filtering, the users need to add flags -low_f in the CBIG_preproc_fslmcflirt_outliers.csh, the function will then perform bandstop/low-pass filtering, depending on if -high_f is used
  - If -low_f is passed in but -high_f is empty, we use low pass filter and -low_f is the stop frequency;
  - If both -low_f and -high_f are provided we use bandstop filter and the stopband is [low_f,high_f]
  - If -low_f is not provided we don't perform respiratory pseudomotion filtering
* Example for bandstop filtering: $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fslmcflirt_outliers.csh 
	-s ABCD_sub1 -d ~/storage/FMRI_preprocess -bld '002 003' -BOLD_stem _rest_skip4_stc -nframe 0 
	-FD_th 0.2 -DVARS 50 -discard-run 50 -rm-seg 5 -low_f 0.31 -high_f 0.43
* Example for lowpass filtering: $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fslmcflirt_outliers.csh 
	-s ABCD_sub1 -d ~/storage/FMRI_preprocess -bld '002 003' -BOLD_stem _rest_skip4_stc -nframe 0 
	-FD_th 0.2 -DVARS 50 -discard-run 50 -rm-seg 5 -low_f 0.1
	
## References

Fleming, S., Thompson, M., Stevens, R., Heneghan, C., Plüddemann, A., Maconochie, I., Tarassenko, L., and Mant, D., 2011. Normal ranges of heart rate and respiratory rate in children from birth to 18 years of age: a systematic review of observational studies. The Lancet, 377 (9770), 1011-1018.


Power, J.D., Lynch, C.J., Silver, B.M., Dubin, M.J., Martin, A., and Jones, R.M., 2019. Distinctions among real and apparent respiratory motions in human fMRI data. NeuroImage, 201, 116041.


Fair, D.A., Miranda-Dominguez, O., Snyder, A.Z., Perrone, A., Earl, E.A., Van, A.N., Koller, J.M., Feczko, E., Tisdall, M.D., van der Kouwe, A., Klein, R.L., Mirro, A.E., Hampton, J.M., Adeyemo, B., Laumann, T.O., Gratton, C., Greene, D.J., Schlaggar, B.L., Hagler, D.J., Jr, Watts, R., Garavan, H., Barch, D.M., Nigg, J.T., Petersen, S.E., Dale, A.M., Feldstein-Ewing, S.W., Nagel, B.J., and Dosenbach, N.U.F., 2020. Correction of respiratory artifacts in MRI head motion estimates. NeuroImage, 208, 116400.


Gratton, C., Dworetsky, A., Coalson, R.S., Adeyemo, B., Laumann, T.O., Wig, G.S., Kong, T.S., Gratton, G., Fabiani, M., Barch, D.M., Tranel, D., Miranda-Dominguez, O., Fair, D.A., Dosenbach, N.U.F., Snyder, A.Z., Perlmutter, J.S., Petersen, S.E., and Campbell, M.C., 2020. Removal of high frequency contamination from motion estimates in single-band fMRI saves data without biasing functional connectivity. NeuroImage, 217, 116866.
