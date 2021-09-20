# Input data used for Parametric Mean Field Model Analysis
* `Desikan_input` this folder contains the input data for pMFM analysis in Desikan parcellation
* `Schaefer100_input` this folder contains the input data for pMFM analysis in Schaefer parcellation
* Those data are used in the folder `part1_pMFM_main` and `part2_pMFM_control_analysis/XXXXX` by creating a symbolic link under each analysis named `input`
* Users can also delete the symbolic links and create a new one for users' own data in the folder `part1_pMFM_main` and `part2_pMFM_control_analysis/XXXXX`

## Data in `Desikan_input`
* `sc_train.csv` Group level SC matrix of training set
* `sc_vali.csv` Group level SC matrix of validation set
* `sc_test.csv` Group level SC matrix of test set
* `subject_list.mat` Subject index of training, validation and test sets

* `fc_train.csv` Averaged empirical FC matrix of training set
* `fc_vali.csv` Averaged empirical FC matrix of validation set
* `fc_test.csv` Averaged empirical FC matrix of test set

* `fcd_train.mat` Averaged empirical FCD CDF of training set
* `fcd_vali.mat` Averaged empirical FCD CDF of validation set
* `fcd_test.mat` Averaged empirical FCD CDF of test set
* `fcd_test_high_window.mat` Averaged empirical FCD CDF using 125 TR sliding window of test set
* `fcd_test_low_window.mat` Averaged empirical FCD CDF using 43 TR sliding window of test set

* `GS_sort.mat` Variance of global signal for all HCP subjects. 
	* This data is copied from [Orban2020_tod](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/preprocessing/Orban2020_tod)

* `rsfc_gradient.csv` Group level first principal gradient of training set
* `rsfc_gradient_pc2.csv` Group level second principal gradient of training set
* `rsfc_gradient_spinned.csv` Spinned group level first principal gradient of training set

* `myelin.csv` Group level T1w/T2w myelin estimation of training set
* `myelin_spinned.csv` Spinned group level T1w/T2w myelin estimation of training set

* `intersubject_variability.mat` Group level inter-subject function connectivity variability of training set

* `structual_gradient1.mat` Group level first gradient decomposition of structural covariance of training set

* `run_label_testset.mat` Indicator of each run belonging to which subject. In detail, each run has a label of subject index. This information will used for bootstrapping.

* `Desikan_68_region_description.mat` Region description of the Desikan 68

* `gene_data.mat`,`gene_spin_data.mat`,`gene_expression_data.csv`,`gene_gradient1.mat` are gene related data
	* The scripts to generate gene expression data are from [Anderson 2020](https://www.nature.com/articles/s41467-020-16710-x)
	* `gene_data.mat` Parcellated PVALB, SST and gene expression first principal component 
	* `gene_spin_data.mat` 1000 spinned data for parcellated PVALB, SST and gene expression first principal component
	* `gene_expression_data.csv` Gene expression data
	* `gene_gradient1.mat` First principal component of gene expression data

* `individual_input` folder contains the 12 individual data for training, validation and test
	
	
## Data in `Schaefer100_input`
* `sc_train.csv` Group level SC matrix of training set
* `sc_vali.csv` Group level SC matrix of validation set
* `sc_test.csv` Group level SC matrix of test set

* `fc_train.csv` Averaged empirical FC matrix of training set
* `fc_vali.csv` Averaged empirical FC matrix of validation set
* `fc_test.csv` Averaged empirical FC matrix of test set

* `fcd_train.mat` Averaged empirical FCD CDF of training set
* `fcd_vali.mat` Averaged empirical FCD CDF of validation set
* `fcd_test.mat` Averaged empirical FCD CDF of test set
* `fcd_test_high_window.mat` Averaged empirical FCD CDF using 125 TR sliding window of test set
* `fcd_test_low_window.mat` Averaged empirical FCD CDF using 43 TR sliding window of test set

* `rsfc_gradient` Group level first principal gradient of training set

* `Schaefer2018_100Parcels_17Networks_order.dlabel.nii` Schaefer 100 parcellation file

* `myelin.csv` Group level T1w/T2w myelin estimation of training set

* `run_label_testset.mat` Indicator of each run belonging to which subject. In detail, each run has a label of subject index. This information will used for bootstrapping.

* `Schaefer_100_region_description.mat` Region description of the Desikan 68

* `gene_data.mat`,`gene_spin_data.mat`,`gene_expression_data.csv` are gene related data
	* The scripts to generate gene expression data are from [Anderson 2020](https://www.nature.com/articles/s41467-020-16710-x)
	* `gene_data.mat` Parcellated PVALB, SST and gene expression first principal component 
	* `gene_spin_data.mat` 1000 spinned data for parcellated PVALB, SST and gene expression first principal component 
	* `gene_expression_data.csv` Gene expression data
	
