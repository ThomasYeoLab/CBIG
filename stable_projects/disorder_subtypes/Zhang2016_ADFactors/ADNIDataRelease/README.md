## Want to extract factor compositions for non-ADNI participants?

Please go to folder `inferNew` for instructions and code.

----

## Looking for the probabilistic atrophy maps?

Please go to folders `inferNew/files_LDA/model_K2`, `inferNew/files_LDA/model_K3` and `inferNew/files_LDA/model_K4` for the two-, three- and four-factor models, respectively. You can also view or download all of these atrophy maps online at [http://neurovault.org/collections/1917/](http://neurovault.org/collections/1917/).

----

## Factor Compositions of ADNI Participants

Spreadsheet `atrophyFactorCompositions_ADNI1_810bl-560m24.csv` includes 810 ADNI1-enrolled subjects’ baseline factor compositions and their factor compositions 24 months after baseline (N = 560). Therefore, it has 1370 rows (excluding the header row).

### Column Definitions
It has eight columns: `RID`, `VISCODE2`, `ImageUID`, `Phase`, `ScanDate`, `Cortical_Prob`, `Temporal_Prob`, and `Subcortical_Prob`.

1. `RID` (roster ID). This column corresponds to the `RID` column of, e.g., the diagnosis file `DXSUM_PDXCONV_ADNIALL.csv` (downloadable from the ADNI website). It also matches with the last four digits of the `SubjectID` column of `MPRAGEMETA.csv` (downloadable from the ADNI website).
2. `VISCODE2` (visit code). This column corresponds to the `VISCODE2` column of the diagnosis file `DXSUM_PDXCONV_ADNIALL.csv`. Possible values in this spreadsheet are `bl`, standing for baseline (or screening), and `m24`, standing for 24 months after baseline.
3. `ImageUID` (image ID). This column corresponds to the `ImageUID` column of `MPRAGEMETA.csv`. Using `ImageUID`, You can uniquely identify an image within the whole ADNI dataset. **Note --**There are seven images on our server that are no longer available on the ADNI website, so their `ImageUID` values are `NA`.
4. `Phase` (study phase). In this spreadsheet, its value is always `ADNI1`.
5. `ScanDate` (scan date). This column corresponds to the `ScanDate` column of `MPRAGEMETA.csv`. Note that this date is different from, but usually very close to, the examination date (the `EXAMDATE` column of the diagnosis file `DXSUM_PDXCONV_ADNIALL.csv`).  **Note --** The seven images mentioned in 3. have `NA` in this column too.
6. `Cortical_Prob` (cortical factor probability).
7. `Temporal_Prob` (temporal factor probability).
8. `Subcortical_Prob` (subcortical factor probability).

### Loading the Spreadsheet into MATLAB for Further Analyses

We provide the following code snippet that loads the spreadsheet into MATLAB. It can be easily translated to another language such as R.

```
% Load
fID = fopen('atrophyFactorCompositions_ADNI1_810bl-560m24.csv');
data = textscan(fID, '%d %s %s %s %s %f %f %f', 'Delimiter', ',', 'HeaderLines', 1);
fclose(fID);

% Extract
RID = data{1};
VISCODE2 = data{2};
ImageUID = data{3};
Phase = data{4};
ScanDate = data{5};
Cortical_Prob = data{6};
Temporal_Prob = data{7};
Subcortical_Prob = data{8};

```

While parsing the spreadsheet, please bear in mind the `NA` values in `ImageUID` and `ScanDate` columns (see the previous section).
