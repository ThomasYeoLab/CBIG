## Reference

Orban C, Kong R, Li J, Chee MWL, Yeo BTT. [Time of day is associated with paradoxical reductions in global signal fluctuation and functional connectivity](https://doi.org/10.1101/653899). In press

## Background

The brain exhibits substantial diurnal variation in physiology and function but neuroscience studies rarely report or consider the effects of time of day. Here, we examined variation in resting-state fMRI in around 900 subjects scanned between 8am to 10pm on two different days. Multiple studies across animals and humans have demonstrated that the brain’s global signal amplitude (henceforth referred to as “fluctuation”) increases with decreased arousal. Thus, in accord with known circadian variation in arousal, we hypothesised that global signal fluctuation would be lowest in the morning, increase in the mid-afternoon and dip in the early evening. Instead, we observed a cumulative decrease (22% between 9am to 9pm) in global signal fluctuation as the day progressed. To put the magnitude of this decrease in context, we note that task-evoked fMRI responses are typically in the order of 1% to 3%. Respiratory variation also decreased with time of day, although control analyses suggested that this did not account for the reduction in GS fluctuation. Finally, time of day was associated with marked decreases in resting state functional connectivity across the whole brain. The magnitude of decrease was significantly stronger than associations between functional connectivity and behaviour (e.g., fluid intelligence). These findings reveal unexpected effects of time of day on the resting human brain, which challenge the prevailing notion that the brain’s global signal reflects mostly arousal and physiological artefacts. We conclude by discussing potential mechanisms for the observed diurnal variation in resting brain activity and the importance of accounting for time of day in future studies.

## Usage

We have added two csv files in the `data_release` folder:
+ `HCP_S1200_physio_data_summary_2020_02_11.csv` 
This csv file contains various run-level summary metrics, including our visual quality assessment (pass/fail) of all respiratory and pulse time-series in the entire S1200 dataset.

+ `HCP_S1200_dictionary_2020_02_11.csv`
This csv file contains a dictionary which defines all variables in `HCP_S1200_physio_data_summary_2020_02_11.csv`

We will soon update this repo with further details.

## Download

To download the data, you can either

+ visit this link: [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.18.1-Update_stable_project_unit_test](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.18.1-Update_stable_project_unit_test)

or

+ run the following command, if you have Git installed

```
git checkout -b Orban2020_tod v0.18.1-Update_stable_project_unit_test
```

## Updates

+ Release v0.18.1 (20/01/2021): Update config file.
+ Release v0.16.2 (11/02/2020): Updated README file for Orban2020_tod and changed labeling convention from 'std' to 'SD' in csv files
+ Release v0.16.0 (19/12/2019): Initial release of Orban2020_tod

## Bugs and Questions

Please contact Csaba Orban at csaba.orban@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
